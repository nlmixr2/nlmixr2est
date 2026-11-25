#include <Rcpp.h>
#include <string>
#include <cstdio>
#include <sys/stat.h>
#include <chrono>
#include "stopFile.h"

static std::string _stopFile;
static bool _stopFileOn = false;
static bool _stopRequested = false;
static bool _stopped = false;

// maxTime= (#1514): a wall-clock limit that requests a stop once, like the file
static bool _timeOn = false;
static bool _timeArmed = false;
static bool _timedOut = false;
static double _maxTime = 0.0;
static std::chrono::steady_clock::time_point _startTime;

bool nmStopFileRequested() {
  return _stopRequested;
}

const char *nmStopMessage() {
  return _timedOut ? "stopped at the time limit" : "stopped early by user request";
}

bool nmStopFilePoll() {
  if (_stopRequested) return true;
  if (_timeArmed) {
    std::chrono::duration<double> el = std::chrono::steady_clock::now() - _startTime;
    if (el.count() >= _maxTime) {
      // one-shot: a later loop of the same run is not stopped again
      _timeArmed = false;
      _stopRequested = _stopped = _timedOut = true;
      return true;
    }
  }
  if (!_stopFileOn) return false;
  struct stat st;
  if (stat(_stopFile.c_str(), &st) != 0) return false;
  // removed so a re-run in the same folder is not stopped at once
  std::remove(_stopFile.c_str());
  // the run is already stopping: a later loop is not stopped by the time limit
  _timeArmed = false;
  _stopRequested = _stopped = true;
  return true;
}

void nmStopFileConsume() {
  _stopRequested = false;
}

// Sets the stop file polled during estimation ("" turns polling off) and clears
// any earlier request; returns the previous file.
//[[Rcpp::export(".nmStopFileSet")]]
std::string nmStopFileSet_(std::string file) {
  std::string old = _stopFileOn ? _stopFile : std::string("");
  _stopFile = file;
  _stopFileOn = !file.empty();
  _stopRequested = _stopped = false;
  return old;
}

// Starts the wall-clock limit of maxTime seconds from now (NA turns it off) and
// clears an earlier time-out.
//[[Rcpp::export(".nmStopTimeSet")]]
void nmStopTimeSet_(double maxTime) {
  _timeOn = _timeArmed = !ISNAN(maxTime);
  _maxTime = _timeOn ? maxTime : 0.0;
  _timedOut = false;
  _startTime = std::chrono::steady_clock::now();
}

// Whether a stop file or a time limit is being polled.
//[[Rcpp::export(".nmStopActive")]]
bool nmStopActive_() {
  return _stopFileOn || _timeOn;
}

// Whether the time limit ended a loop since the last .nmStopTimeSet().
//[[Rcpp::export(".nmStopTimedOut")]]
bool nmStopTimedOut_() {
  return _timedOut;
}

// The message of the stop that was requested.
//[[Rcpp::export(".nmStopMessage")]]
std::string nmStopMessage_() {
  return std::string(nmStopMessage());
}

// The stop file being polled ("" when none).
//[[Rcpp::export(".nmStopFileGet")]]
std::string nmStopFileGet_() {
  return _stopFileOn ? _stopFile : std::string("");
}

// Polls the stop file; TRUE once a stop has been requested.
//[[Rcpp::export(".nmStopFilePoll")]]
bool nmStopFilePoll_() {
  return nmStopFilePoll();
}

// Whether a stop has been requested, without polling.
//[[Rcpp::export(".nmStopFileRequested")]]
bool nmStopFileRequested_() {
  return nmStopFileRequested();
}

// Whether a stop request ended a loop since the last .nmStopFileSet().
//[[Rcpp::export(".nmStopFileStopped")]]
bool nmStopFileStopped_() {
  return _stopped;
}

// Marks a stop request as handled by the loop it ended.
//[[Rcpp::export(".nmStopFileConsume")]]
void nmStopFileConsume_() {
  nmStopFileConsume();
}
