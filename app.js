import { makeApi, parseRecovery, validateFiles } from "./api.js";

const $ = selector => document.querySelector(selector);
const show = (selector, visible = true) => { $(selector).hidden = !visible; };
const text = (selector, value) => { $(selector).textContent = value; };
const number = value => Number(value ?? 0).toLocaleString();
const fileSize = bytes => { const power = bytes >= 1024**3 ? 3 : bytes >= 1024**2 ? 2 : 1; return `${(bytes / 1024**power).toLocaleString(undefined, { maximumFractionDigits: 2 })} ${["", "KiB", "MiB", "GiB"][power]}`; };
const terminal = new Set(["completed", "failed", "cancelled"]);
let api, credential = parseRecovery(location.hash), currentJob = null, result = null;
let files = [], busy = false, ready = false, example = false, poll = null, generation = 0;
let maxBytes = 16 * 1024**3, pageRequest = 0, searchTimer;
function notice(message = "") { text("#notice", message); show("#notice", Boolean(message)); }
function friendly(error) { return error.name === "TimeoutError" || error.name === "AbortError" ? "The service did not respond in time. Use Refresh status to reconnect; your job may still be running." : error instanceof TypeError ? "The service cannot be reached. Check your connection and try again; your job may still be running." : error.message; }
function mode() { return $("input[name=mode]:checked").value; }
function enableSubmit() {
  document.querySelectorAll("#analysis-form input, #analysis-form select, #analysis-form button").forEach(control => { control.disabled = busy; });
  $("#low-abundance").disabled = busy || mode() !== "reads";
  $("#submit").disabled = busy || !ready || Boolean(credential);
}
function options() {
  const values = { boundary_mode: $("#boundary-mode").value, low_abundance: mode() === "reads" && $("#low-abundance").checked };
  for (const kind of ["repeat", "spacer"]) {
    const lo = Number($(`#${kind}-min`).value), hi = Number($(`#${kind}-max`).value);
    if (![lo, hi].every(value => Number.isInteger(value) && value >= 0 && value <= 500) || (lo && hi && lo > hi)) throw new Error(`Check the ${kind} length range. Minimum must not exceed maximum.`);
    values[`${kind}_min_length`] = lo; values[`${kind}_max_length`] = hi;
  }
  return values;
}
function fileSummary() {
  text("#file-summary", files.length ? files.map(file => `${file.name} (${fileSize(file.size)})`).join("; ") : "");
  show("#file-summary", files.length > 0);
}
function changeMode() {
  const reads = mode() === "reads";
  $("#files").multiple = reads;
  $("#files").accept = reads ? ".fastq,.fq,.fasta,.fa,.fna,.gz" : ".txt,.gz";
  $("#low-abundance").disabled = !reads;
  text("#input-label", reads ? "Upload reads" : "Upload retained cycles");
  text("#input-hint", reads ? "FASTA / FASTQ · 1–2 files · .gz accepted" : "cycles.txt · 1 file · .gz accepted");
  files = []; $("#files").value = ""; fileSummary();
}
async function service() {
  if (!api) return;
  try {
    const [health, config] = await Promise.all([api.request("/health"), api.request("/config")]);
    ready = health.status === "ready";
    maxBytes = config.max_input_bytes;
    text("#limits", `${maxBytes >= 1024**3 ? maxBytes / 1024**3 + " GiB" : maxBytes / 1024**2 + " MiB"} max · ${config.retention_hours / 24}-day retention`);
  } catch { ready = false; }
  $("#service").classList.toggle("online", ready);
  text("#service span:last-child", ready ? "Service ready" : "Service unavailable");
  show("#check-service", !ready);
  enableSubmit();
}
function recovery() {
  if (!credential) return;
  const link = new URL(location.href); link.search = ""; link.hash = `job=${credential.id}.${credential.token}`;
  history.replaceState(null, "", link);
  $("#recovery-link").value = link.href;
  text("#job-id", `Job ${credential.id.slice(0, 8)}`);
}
function displayJob(job) {
  currentJob = job; example = false;
  show("#analysis-view", false); show("#job-view");
  text("#job-state", job.status[0].toUpperCase() + job.status.slice(1));
  $("#job-state").dataset.state = job.status;
  text("#job-message", job.message || "Uploading input files.");
  if (job.expires) text("#expires", `Available until ${new Date(job.expires * 1000).toLocaleString()}`);
  show("#cancel", !terminal.has(job.status));
  $("#cancel").disabled = busy || job.status === "cancelling";
  const order = ["uploading", "queued", "running", "completed"];
  document.querySelectorAll(".job-steps li").forEach(item => {
    item.classList.toggle("done", order.indexOf(item.dataset.step) < order.indexOf(job.status));
    item.classList.toggle("active", item.dataset.step === job.status);
  });
  recovery();
}
async function refresh() {
  if (!credential || busy) return;
  clearTimeout(poll);
  const session = credential, version = generation;
  try {
    const job = await api.request(`/jobs/${session.id}`, { credential: session });
    if (version !== generation) return;
    notice(); displayJob(job);
    if (job.status === "completed" && !result) {
      const output = await api.request(`/jobs/${session.id}/result`, { credential: session });
      if (version !== generation) return;
      renderResult(output, false);
    }
    if (!terminal.has(job.status)) poll = setTimeout(refresh, 3000);
  } catch (error) {
    if (version !== generation) return;
    notice(friendly(error));
    if (![404, 410].includes(error.status)) poll = setTimeout(refresh, 15000);
  }
}
function reset() {
  if (busy) return;
  if (credential && !confirm("Have you saved your recovery link? Starting another analysis removes this link from the address bar. An active job will continue running.")) return;
  generation++; pageRequest++; clearTimeout(searchTimer); clearTimeout(poll); credential = null; currentJob = null; result = null; example = false;
  history.replaceState(null, "", location.pathname + location.search);
  show("#job-view", false); show("#result-view", false); show("#analysis-view"); notice(); enableSubmit();
  $("#analysis-heading").setAttribute("tabindex", "-1"); $("#analysis-heading").focus({ preventScroll: true });
  $("#analysis-heading").scrollIntoView({ block: "start" });
}
function element(tag, value, className) {
  const node = document.createElement(tag); if (value !== undefined) node.textContent = value; if (className) node.className = className; return node;
}
function renderGroups() {
  const container = $("#groups"); container.replaceChildren();
  const query = $("#group-search").value.trim().toUpperCase();
  const rows = new Map();
  for (const row of result.entries) { if (!rows.has(row.array_id)) rows.set(row.array_id, []); rows.get(row.array_id).push(row); }
  let found = 0;
  for (const group of result.arrays) {
    const entries = rows.get(group.array_id) || [];
    if ((example || !result.pagination) && query && ![group.array_id, group.repeat_consensus, ...entries.map(row => row.spacer)].some(value => String(value || "").toUpperCase().includes(query))) continue;
    found++;
    const details = element("details", undefined, "group");
    details.open = result.arrays.length === 1;
    const summary = element("summary");
    summary.append(element("strong", group.array_id), element("span", `${number(group.entry_count)} entries`, "entry-count"), element("code", group.repeat_consensus, "repeat-sequence"));
    details.append(summary);
    const wrap = element("div", undefined, "table-wrap");
    const table = element("table");
    const caption = element("caption", "Observed repeat–spacer pairs. Row order does not indicate genomic order."); table.append(caption);
    const head = element("thead"), tr = element("tr");
    for (const name of ["Entry", "Repeat variant", "Spacer", "Spacer length"]) { const th = element("th", name); th.scope = "col"; tr.append(th); }
    head.append(tr); table.append(head);
    const body = element("tbody");
    const fill = values => {
      body.replaceChildren();
      for (const row of values) {
        const line = element("tr");
        for (const [i, value] of [row.entry_index, row.repeat_variant, row.spacer, row.spacer?.length].entries()) {
          const cell = element("td"); cell.append(element(i === 1 || i === 2 ? "code" : "span", value ?? "Not available")); line.append(cell);
        }
        body.append(line);
      }
    };
    fill(entries);
    if (!example && group.entries_truncated) {
      let offset = 0;
      const controls = element("div",undefined,"pagination"), prev = element("button","Previous entries","secondary"), next = element("button","Next entries","secondary"), label = element("span");
      const update = () => { prev.disabled = offset === 0; next.disabled = offset + 100 >= Number(group.entry_count); label.textContent = `${number(offset+1)}–${number(Math.min(offset+100,Number(group.entry_count)))} of ${number(group.entry_count)}`; };
      const load = async target => {
        prev.disabled = next.disabled = true;
        try {
          const data = await api.request(`/jobs/${credential.id}/groups/${encodeURIComponent(group.array_id)}/entries?offset=${target}`,{credential});
          offset = data.offset; fill(data.entries); notice();
        } catch(error) { notice(friendly(error)); }
        finally { update(); }
      };
      prev.addEventListener("click",()=>load(Math.max(0,offset-100))); next.addEventListener("click",()=>load(offset+100));
      update(); controls.append(prev,label,next); details.append(controls);
    }
    table.append(body); wrap.append(table); details.append(wrap); container.append(details);
  }
  const paging = result.pagination;
  show("#group-pages", !example && Boolean(paging));
  if (paging) {
    text("#group-page-status", `${number(paging.total ? paging.offset+1 : 0)}–${number(Math.min(paging.offset+paging.limit,paging.total))} of ${number(paging.total)} groups`);
    $("#groups-prev").disabled = paging.offset === 0;
    $("#groups-next").disabled = paging.offset + paging.limit >= paging.total;
  }
  if (!found) container.append(element("p", query ? "No groups match your search." : "No repeat–spacer groups were detected with these settings. This does not establish that CRISPR sequences are absent.", "empty-result"));
}
function renderResult(output, isExample) {
  result = output; example = isExample;
  show("#result-view"); show("#download", !isExample); show("#example-download", isExample);
  text("#result-context", isExample ? "Precomputed example · synthetic cycle dataset" : `Analysis completed in ${number(output.elapsed_seconds)} seconds.`);
  const stats = $("#statistics"); stats.replaceChildren();
  for (const [label, value] of [["Primary groups", output.counts.primary_arrays], ["Repeat–spacer entries", output.counts.primary_entries], ["Read-linked groups", output.counts.read_linked_groups]]) {
    const item = element("div"); item.append(element("strong", value == null ? "—" : number(value)), element("span", label)); stats.append(item);
  }
  text("#provenance", JSON.stringify({ settings: output.settings, source_sha256: output.source_sha256, interpretation: output.interpretation, scoring: output.scoring }, null, 2));
  show("#truncated", output.truncated); $("#group-search").value = ""; renderGroups();
}

$("#files").addEventListener("change", event => { files = [...event.target.files]; fileSummary(); });
document.querySelectorAll("input[name=mode]").forEach(input => input.addEventListener("change", changeMode));
for (const eventName of ["dragover", "drop"]) $("#drop-zone").addEventListener(eventName, event => { event.preventDefault(); if (eventName === "drop" && !busy) { files = [...event.dataTransfer.files]; fileSummary(); } });
$("#use-example").addEventListener("click", async () => {
  try {
    const response = await fetch("examples/cycles.txt"); if (!response.ok) throw new Error("Example could not be loaded.");
    $("input[value=cycles]").checked = true; changeMode();
    files = [new File([await response.blob()], "cycles.txt", { type: "text/plain" })]; fileSummary(); notice();
  } catch (error) { notice(friendly(error)); }
});
$("#view-example").addEventListener("click", async () => {
  try {
    const response = await fetch("examples/result.json"); if (!response.ok) throw new Error("Example could not be loaded.");
    renderResult(await response.json(), true); notice(); $("#result-view").scrollIntoView({ block: "start" });
  } catch (error) { notice(friendly(error)); }
});
$("#analysis-form").addEventListener("submit", async event => {
  event.preventDefault(); if (busy || credential) return;
  try {
    validateFiles(files, mode(), maxBytes); const chosenOptions = options();
    busy = true; enableSubmit(); notice();
    text("#submit", "Submitting…"); show("#upload-status"); show("#upload-progress"); text("#upload-status", "Reserving your job…");
    credential = await api.request("/jobs", { method: "POST", body: { mode: mode(), options: chosenOptions } });
    generation++; result = null; show("#result-view", false); recovery();
    let uploaded = 0; const total = files.reduce((sum, file) => sum + file.size, 0);
    for (let index = 0; index < files.length; index++) {
      text("#upload-status", `Uploading file ${index + 1} of ${files.length}…`);
      await api.upload(credential, index, files[index], amount => { $("#upload-progress").value = Math.min(100, (uploaded + amount) / total * 100); });
      uploaded += files[index].size;
    }
    const job = await api.request(`/jobs/${credential.id}/submit`, { credential, method: "POST" });
    displayJob(job); $("#job-view").scrollIntoView({ block: "start" });
  } catch (error) {
    notice(friendly(error));
    if (credential) displayJob({ id: credential.id, status: "uploading", message: "Submission was interrupted. Refresh to check the job, or cancel it before trying again." });
  } finally {
    busy = false; text("#submit", "Run analysis"); show("#upload-status", false); show("#upload-progress", false); enableSubmit();
    $("#cancel").disabled = currentJob?.status === "cancelling";
    if (credential) poll = setTimeout(refresh, 1000);
  }
});
$("#cancel").addEventListener("click", async () => {
  if (!credential || busy) return;
  $("#cancel").disabled = true;
  try { displayJob(await api.request(`/jobs/${credential.id}/cancel`, { credential, method: "POST" })); notice(); await refresh(); }
  catch (error) { notice(`Cancellation is not confirmed. ${friendly(error)}`); $("#cancel").disabled = false; }
});
$("#copy-link").addEventListener("click", async () => {
  try { await navigator.clipboard.writeText($("#recovery-link").value); text("#copy-link", "Copied"); }
  catch { $("#recovery-link").select(); notice("Copy the selected recovery link using your browser’s copy command."); }
});
$("#download").addEventListener("click", async () => {
  if (!credential || example) return;
  $("#download").disabled = true;
  try {
    const ticket = await api.request(`/jobs/${credential.id}/download-ticket`, {credential,method:"POST"});
    const anchor = element("a"); anchor.href = api.downloadUrl(ticket.path); anchor.rel = "noreferrer"; anchor.click();

  } catch (error) { notice(friendly(error)); }
  finally { $("#download").disabled = false; }
});
async function loadPage(offset = 0) {
  if (example || !credential || !result?.pagination) { renderGroups(); return; }
  const session = credential, version = generation, request = ++pageRequest;
  const query = $("#group-search").value.trim();
  $("#groups-prev").disabled = $("#groups-next").disabled = true;
  try {
    const output = await api.request(`/jobs/${session.id}/result?offset=${offset}&q=${encodeURIComponent(query)}`,{credential:session,timeout:60000});
    if (version !== generation || request !== pageRequest) return;
    result = output; renderGroups(); notice();
  } catch(error) {
    if (version === generation && request === pageRequest) { notice(friendly(error)); renderGroups(); }
  }
}
$("#group-search").addEventListener("input", () => { clearTimeout(searchTimer); pageRequest++; searchTimer = setTimeout(()=>loadPage(),250); });
$("#groups-prev").addEventListener("click",()=>loadPage(Math.max(0,result.pagination.offset-result.pagination.limit)));
$("#groups-next").addEventListener("click",()=>loadPage(result.pagination.offset+result.pagination.limit));
$("#refresh").addEventListener("click", refresh);
$("#check-service").addEventListener("click", service);
$("#new-analysis").addEventListener("click", reset);
$("#nav-analysis").addEventListener("click", () => { if (credential) $("#job-view").scrollIntoView({ block: "start" }); else reset(); });
function theme(dark) { document.documentElement.dataset.theme = dark ? "dark" : "light"; text("#theme", dark ? "Light" : "Dark"); $("#theme").setAttribute("aria-label", `Switch to ${dark ? "light" : "dark"} theme`); }
try { theme(localStorage.getItem("zuck-theme") === "dark"); } catch {}
$("#theme").addEventListener("click", () => { const dark = document.documentElement.dataset.theme !== "dark"; theme(dark); try { localStorage.setItem("zuck-theme", dark ? "dark" : "light"); } catch {} });
window.addEventListener("online", () => { service(); refresh(); });
window.addEventListener("focus", () => { if (credential && !busy) refresh(); });
window.addEventListener("hashchange", () => {
  const next = parseRecovery(location.hash);
  if (next && (next.id !== credential?.id || next.token !== credential?.token)) {
    generation++; credential = next; result = null; show("#result-view", false); displayJob({ status: "unknown", message: "Recovering your job…" }); refresh();
  }
});
try { api = makeApi(window.ZUCK_CONFIG?.apiBase); }
catch (error) { notice(error.message); text("#service span:last-child", "Configuration needed"); }
service(); setInterval(service, 60000);
if (credential) { displayJob({ status: "unknown", message: "Recovering your job…" }); refresh(); }
else if (location.hash.startsWith("#job=")) notice("This recovery link is invalid. Use the complete link saved when you submitted the job.");

$('nav a[href="#help"]').addEventListener("click", () => { $("#help").open = true; });
