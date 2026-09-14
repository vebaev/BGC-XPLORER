'use strict';

function isNetworkFetchError(error) {
  const detail = error && error.message ? error.message : '';
  return detail === 'Load failed' || detail === 'Failed to fetch';
}

function fetchWithNetworkRetry(endpoint, options, attempts, delayMs, fetchImpl, waitImpl) {
  const request = fetchImpl || fetch;
  const wait = waitImpl || function (milliseconds) {
    return new Promise(function (resolve) { window.setTimeout(resolve, milliseconds); });
  };
  return request(endpoint, options).catch(function (error) {
    if (attempts <= 1 || !isNetworkFetchError(error)) {
      throw error;
    }
    return wait(delayMs).then(function () {
      return fetchWithNetworkRetry(endpoint, options, attempts - 1, delayMs, request, wait);
    });
  });
}

if (typeof module !== 'undefined' && module.exports) {
  module.exports = {fetchWithNetworkRetry, isNetworkFetchError};
}
