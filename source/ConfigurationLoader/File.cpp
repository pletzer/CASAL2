/**
 * @file File.cpp
 * @author  Scott Rasmussen (scott.rasmussen@zaita.com)
 * @version 1.0
 * @date 16/11/2012
 * @section LICENSE
 *
 * Copyright NIWA Science �2012 - www.niwa.co.nz
 *
 * $Date: 2008-03-04 16:33:32 +1300 (Tue, 04 Mar 2008) $
 */

// Headers
#include "File.h"

#include <boost/algorithm/string/replace.hpp>
#include <boost/algorithm/string/trim_all.hpp>
#include <boost/algorithm/string/split.hpp>

#include "Loader.h"
#include "Utilities/Logging/Logging.h"
#include "Utilities/String.h"
#include "Utilities/To.h"

// Namespaces
namespace isam {
namespace configuration {

namespace util = isam::utilities;

/**
 * Default Constructor
 */
File::File(Loader* loader)
  : loader_(loader) {
}

/**
 * Destructor
 */
File::~File() {
  if (file_.is_open())
    file_.close();
}

/**
 * Attempt to open our configuration file
 *
 * @param file_name The name of the file to open
 * @return true on success, false on failure
 */
bool File::OpenFile(string file_name) {
  LOG_TRACE();

  file_name_ = file_name;

  file_.open(file_name_.c_str());
  if (file_.fail() || !file_.is_open())
    return false;

  return true;
}

/**
 * Parse the configuration file. Creating objects and loading
 * the parameter objects
 */
void File::Parse() {
  LOG_TRACE();

  if (file_.fail() || !file_.is_open())
    LOG_CODE_ERROR("Unable to parse the configuration file because a previous error has not been reported.\nFile: " << file_name_);

  /**
   * Iterate through our file parsing the contents
   */
  bool      multi_line_comment  = false;
  string    current_line        = "";
  while (getline(file_, current_line)) {
    ++line_number_;

    if (current_line.length() == 0)
      continue;

    /**
     * If we're in a multi-line comment we need to look for the end, or discard the line
     */
    if (multi_line_comment) {
      size_t pos = current_line.find_first_of(CONFIG_MULTI_COMMENT_END);
      if (pos == string::npos)
        continue;

      multi_line_comment = false;
      current_line = current_line.substr(pos + 1, current_line.length() - pos);
    }

    /**
     * Check if we're entering a multi-line comment. Strip any comment parts of the line
     */
    size_t pos = current_line.find_first_of(CONFIG_MULTI_COMMENT_START);
    if (pos != string::npos) {
      multi_line_comment = true;
      current_line = current_line.substr(0, pos);
    }

    /**
     * Check and remove any single-line (end of line) comments
     * e.g <data> #comment
     */
    pos = current_line.find_first_of(CONFIG_SINGLE_COMMENT);
    if (pos != string::npos)
      current_line = current_line.substr(0, pos);

    if (current_line.length() == 0)
      continue;
    LOG_INFO("current_line == " << current_line);

    /**
     * Change tabs to spaces, remove any leading/trailing or multiple spaces
     * so we can be sure the input is nicely formatted
     */
    boost::replace_all(current_line, "\t", " ");
    boost::trim_all(current_line);

    /**
     * Now we need to check if this line is an include line for a new
     * file.
     */
    if (current_line.length() > strlen(CONFIG_INCLUDE)) {
      string possible_include = current_line.substr(0, strlen(CONFIG_INCLUDE));
      possible_include = util::ToLowercase(possible_include);

      if (possible_include == CONFIG_INCLUDE) {
        size_t first_quote = current_line.find_first_of("\"");
        size_t last_quote  = current_line.find_last_of("\"");

        if (first_quote == string::npos || last_quote == string::npos || first_quote == last_quote)
          LOG_ERROR("At line " << line_number_ << " of " << file_name_
              << ": File name arguement to " << CONFIG_INCLUDE << " must be surrounded by quotes");

        string include_name = current_line.substr(first_quote);
        boost::replace_all(include_name, "\"", "");
        LOG_INFO("Loading new configuration file via include " << include_name);

        FilePtr include_file = FilePtr(new File(loader_));
        if (!include_file->OpenFile(include_name))
          LOG_ERROR("At line: " << line_number_ << " of " << file_name_
              << ": Include file '" << include_name << "' could not be opened. Does this file exist?");

        include_file->Parse();
        continue;
      }
    }


    /**
     * At this point everything is standard. We have a simple line of text that we now need to parse. All
     * comments etc have been removed and we've gone through any include_file directives
     */
    FileLine current_file_line;
    current_file_line.file_name_    = file_name_;
    current_file_line.line_number_  = line_number_;
    current_file_line.line_         = current_line;

    loader_->AddFileLine(current_file_line);
  } // while(get_line())
}

} /* namespace configuration */
} /* namespace isam */


















