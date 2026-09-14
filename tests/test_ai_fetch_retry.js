'use strict';

const assert = require('assert');
const {fetchWithNetworkRetry} = require('../scripts/ai_fetch_retry.js');

async function main() {
  let attempts = 0;
  const fetchImpl = async () => {
    attempts += 1;
    if (attempts < 3) {
      throw new Error('Load failed');
    }
    return {ok: true};
  };
  const response = await fetchWithNetworkRetry(
    '/analyze_cluster', {}, 3, 0, fetchImpl, () => Promise.resolve()
  );
  assert.strictEqual(response.ok, true);
  assert.strictEqual(attempts, 3);

  let applicationAttempts = 0;
  await assert.rejects(
    fetchWithNetworkRetry(
      '/analyze_cluster', {}, 3, 0,
      async () => {
        applicationAttempts += 1;
        throw new Error('NVIDIA API error 503');
      },
      () => Promise.resolve(),
    ),
    /NVIDIA API error 503/,
  );
  assert.strictEqual(applicationAttempts, 1);
}

main().catch(error => {
  console.error(error);
  process.exit(1);
});
