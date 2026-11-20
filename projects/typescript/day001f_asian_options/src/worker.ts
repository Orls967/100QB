import { parentPort } from "node:worker_threads";
import { type ChunkTask, type OptionParams, runChunk } from "./core.js";

export interface WorkerRequest {
  task: ChunkTask;
  params: OptionParams;
  seed: number;
  hRel: number;
  chunks: Array<{ idx: number; len: number }>;
}

export interface WorkerResponse {
  results: Array<{ idx: number; acc: number[] }>;
}

if (parentPort) {
  parentPort.on("message", (msg: WorkerRequest) => {
    const results = msg.chunks.map((c) => ({
      idx: c.idx,
      acc: Array.from(runChunk(msg.task, msg.params, c.idx, c.len, msg.seed, msg.hRel)),
    }));
    const response: WorkerResponse = { results };
    parentPort?.postMessage(response);
  });
}
