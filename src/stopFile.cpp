#include <Rcpp.h>
#include <string>
#include <cstdio>
#include <sys/stat.h>
#include "stopFile.h"

static std::string _stopFile;
static bool _stopFileOn = false;
static bool _stopRequested = false;
static bool _stopped = false;

bool nmStopFileRequested() {
  return _stopRequested;
}

bool nmStopFilePoll() {
  if (_stopRequested) return true;
  if (!_stopFileOn) return false;
  struct stat st;
  if (stat(_stopFile.c_str(), &st) != 0) return false;
  // removed so a re-run in the same folder is not stopped at once
  std::remove(_stopFile.c_str());
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
