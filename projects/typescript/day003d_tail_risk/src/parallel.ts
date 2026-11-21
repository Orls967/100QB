import { existsSync } from "node:fs";
import os from "node:os";
import { Worker } from "node:worker_threads";
import { executeChunk, type ChunkOutput, type ParallelTaskSpec } from "./core_tasks.js";
import type { WorkerRequest, WorkerResponse } from "./worker.js";

function resolveWorkerUrl(): URL {
  const sameDirJs = new URL("./worker.js", import.meta.url);
  if (existsSync(sameDirJs)) {
    return sameDirJs;
  }
  return new URL("../dist/worker.js", import.meta.url);
}

/**
 * Execute `totalItems` deterministically across fixed chunks of size `chunkSize`.
 * Results are collected and returned strictly ordered by chunkIdx (0..nChunks-1),
 * guaranteeing bit-identical outputs for workers = 1 vs workers = N.
 */
export async function runParallelTask(
  task: ParallelTaskSpec,
  totalItems: number,
  chunkSize: number,
  mainSeed: number,
  salt: number,
  workers?: number,
): Promise<{ outputs: ChunkOutput[]; workerCount: number }> {
  const nChunks = Math.ceil(totalItems / chunkSize);
  const workerCount = Math.min(
    workers && workers > 0 ? workers : os.availableParallelism(),
    Math.max(nChunks, 1),
  );

  if (workerCount <= 1) {
    const outputs: ChunkOutput[] = new Array(nChunks);
    for (let c = 0; c < nChunks; c += 1) {
      const start = c * chunkSize;
      const len = Math.min(totalItems - start, chunkSize);
      outputs[c] = executeChunk(task, mainSeed, salt, c, len);
    }
    return { outputs, workerCount: 1 };
  }

  const buckets: Array<Array<{ idx: number; len: number }>> = Array.from(
    { length: workerCount },
    () => [],
  );
  for (let c = 0; c < nChunks; c += 1) {
    const start = c * chunkSize;
    const len = Math.min(totalItems - start, chunkSize);
    buckets[c % workerCount]?.push({ idx: c, len });
  }

  const workerUrl = resolveWorkerUrl();
  const outputs: ChunkOutput[] = new Array(nChunks);

  await Promise.all(
    buckets.map(
      (chunks) =>
        new Promise<void>((resolve, reject) => {
          const worker = new Worker(workerUrl);
          worker.once("message", (resp: WorkerResponse) => {
            for (const item of resp.results) {
              outputs[item.idx] = item.output;
            }
            void worker.terminate();
            resolve();
          });
          worker.once("error", (err) => {
            void worker.terminate();
            reject(err);
          });
          const req: WorkerRequest = { task, mainSeed, salt, chunks };
          worker.postMessage(req);
        }),
    ),
  );

  return { outputs, workerCount };
}
