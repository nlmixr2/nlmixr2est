#ifndef NLMIXR2EST_STOP_FILE_H
#define NLMIXR2EST_STOP_FILE_H

#include <exception>

// A graceful early stop requested by creating the control's `stopFile` (#1507).
// nmStopFilePoll() checks for the file (removing it when found) and returns
// whether a stop has been requested.  The request holds until the loop it ends
// calls nmStopFileConsume(), so a later fit of the same run is not stopped.
// A control's `maxTime` requests the stop the same way, once its time is up
// (#1514); nmStopMessage() says which stop was requested.
// Call these only from the main thread.
bool nmStopFilePoll();
const char *nmStopMessage();
bool nmStopFileRequested();
void nmStopFileConsume();

// Thrown out of an optimizer's objective callback to end the search.
struct nmStopRequest : public std::exception {
  const char *what() const noexcept override { return "stopped early by user request"; }
};

#endif
