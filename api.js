export function parseRecovery(hash) {
  const match = /^#job=([0-9a-f]{32})\.([A-Za-z0-9_-]{43})$/.exec(hash);
  return match ? { id: match[1], token: match[2] } : null;
}

export function validateFiles(files, mode, maxBytes) {
  if (!files.length) throw new Error("Choose an input file first.");
  if (files.length > (mode === "reads" ? 2 : 1)) throw new Error(mode === "reads" ? "Upload one or two read files." : "Upload one cycles.txt file.");
  if (files.some(file => !file.size)) throw new Error("Input files must not be empty.");
  if (files.reduce((sum, file) => sum + file.size, 0) > maxBytes) throw new Error(`Combined inputs exceed the ${maxBytes >= 1024**3 ? maxBytes / 1024**3 + " GiB" : maxBytes / 1024**2 + " MiB"} upload limit.`);
}

export function makeApi(base) {
  const parsed = new URL(base);
  if (parsed.username || parsed.password || parsed.search || parsed.hash || (parsed.protocol !== "https:" && !(parsed.protocol === "http:" && ["127.0.0.1", "localhost"].includes(parsed.hostname)))) throw new Error("The API address is not configured correctly.");
  base = base.replace(/\/$/, "");
  async function request(path, { credential, method = "GET", body, binary = false, timeout = 15000 } = {}) {
    const headers = {};
    if (credential) headers.Authorization = `Bearer ${credential.token}`;
    if (body !== undefined) headers["Content-Type"] = "application/json";
    const response = await fetch(base + "/api/v1" + path, { method, body: body === undefined ? undefined : JSON.stringify(body), headers, credentials: "omit", referrerPolicy: "no-referrer", cache: "no-store", signal: AbortSignal.timeout(timeout) });
    if (!response.ok) {
      const error = await response.json().catch(() => ({}));
      const failure = new Error(typeof error.detail === "string" ? error.detail : `Service request failed (${response.status}). Try again.`);
      failure.status = response.status;
      throw failure;
    }
    return binary ? response.blob() : response.json();
  }
  function upload(credential, slot, file, onProgress) {
    return new Promise((resolve, reject) => {
      const xhr = new XMLHttpRequest();
      xhr.open("PUT", `${base}/api/v1/jobs/${credential.id}/input/${slot}`);
      xhr.setRequestHeader("Authorization", `Bearer ${credential.token}`);
      xhr.setRequestHeader("Content-Type", "application/octet-stream");
      xhr.timeout = 14400000;
      xhr.upload.onprogress = event => onProgress(event.loaded);
      xhr.onerror = () => reject(new Error("Upload connection lost. Your recovery link is saved below; cancel this job before trying again."));
      xhr.ontimeout = () => reject(new Error("Upload timed out. Cancel this job and try a smaller input."));
      xhr.onload = () => {
        let body = {};
        try { body = JSON.parse(xhr.responseText); } catch {}
        if (xhr.status >= 200 && xhr.status < 300) resolve(body);
        else reject(new Error(typeof body.detail === "string" ? body.detail : "Upload failed. Cancel this job and try again."));
      };
      xhr.send(file);
    });
  }
  return { request, upload, downloadUrl: path => base + "/api/v1" + path };
}
