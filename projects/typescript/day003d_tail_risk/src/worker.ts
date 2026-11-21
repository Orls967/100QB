import { parentPort, type TransferListItem } from "node:worker_threads";
import { executeChunk, type ChunkOutput, type ParallelTaskSpec } from "./core_tasks.js";

export interface WorkerRequest {
  task: ParallelTaskSpec;
  mainSeed: number;
  salt: number;
  chunks: Array<{ idx: number; len: number }>;
}

export interface WorkerResponse {
  results: Array<{ idx: number; output: ChunkOutput }>;
}

if (parentPort) {
  parentPort.on("message", (msg: WorkerRequest) => {
    const results: Array<{ idx: number; output: ChunkOutput }> = [];
    const transfers: TransferListItem[] = [];

    for (const c of msg.chunks) {
      const output = executeChunk(msg.task, msg.mainSeed, msg.salt, c.idx, c.len);
      if (output.kind === "pbh_sample" && output.samples.buffer instanceof ArrayBuffer) {
        transfers.push(output.samples.buffer);
      } else if (output.kind === "z2_mc" && output.z2Flat.buffer instanceof ArrayBuffer) {
        transfers.push(output.z2Flat.buffer);
      }
      results.push({ idx: c.idx, output });
    }

    const response: WorkerResponse = { results };
    parentPort?.postMessage(response, transfers);
  });
}
