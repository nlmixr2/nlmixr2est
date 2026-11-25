#ifndef NLMIXR2EST_STOP_FILE_H
#define NLMIXR2EST_STOP_FILE_H

#include <exception>

// A graceful early stop requested by creating the control's `stopFile` (#1507).
// nmStopFilePoll() checks for the file (removing it when found) and returns
// whether a stop has been requested; the request is sticky until the next
// nmStopFileSet().  Call it only from the main thread.
bool nmStopFilePoll();
bool nmStopFileRequested();

// Thrown out of an optimizer's objective callback to end the search.
struct nmStopRequest : public std::exception {
  const char *what() const noexcept override { return "stopped early by user request"; }
};

#endif
