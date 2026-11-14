// Sequence evidence is shared by consensuses, repeat variants and every spacer.
export function evidenceUrl(value) {
  try {
    const url = new URL(value);
    return url.protocol === "https:" && url.hostname === "www.ncbi.nlm.nih.gov" &&
      /^\/(nuccore|assembly)\//.test(url.pathname) && !url.username && !url.password ? url.href : null;
  } catch { return null; }
}

export function evidenceLabel(match, kind) {
  const accession = match.accession || match.subject_accession || "Source record";
  if (kind === "protospacer") {
    return `${accession}:${match.subject_start}–${match.subject_end} (${match.relative_strand}) · ${match.substitutions} ${match.substitutions === 1 ? "substitution" : "substitutions"} · ${Math.round(match.query_coverage * 100)}% coverage`;
  }
  const coordinates = match.start != null && match.end != null ? `:${match.start}–${match.end}` : "";
  const sourceKind = {consensus:"repeat consensus", repeat:"repeat", spacer:"spacer", array_projection:"source array"}[match.source_kind] || "sequence";
  return `${accession}${coordinates} · ${match.category || "Catalogue"} · ${sourceKind} (${match.relative_strand})`;
}

const node = (tag, text, cls) => {
  const result = document.createElement(tag);
  if (text !== undefined) result.textContent = text;
  if (cls) result.className = cls;
  return result;
};

export function sequenceEvidence(sequence, annotation, loadMore, onError) {
  const wrap = node("div", undefined, "sequence-evidence");
  const code = node("code", sequence);
  if (annotation?.database_count > 0) code.classList.add("database-match");
  if (annotation?.protospacer_count > 0) code.classList.add("protospacer-match");
  wrap.append(code);
  if (!annotation) return wrap;
  for (const kind of ["database", "protospacer"]) {
    const count = Number(annotation[kind + "_count"] || 0);
    const label = kind === "database" ? "catalogue" : "candidate protospacer";
    if (!count) {
      if (kind === "protospacer" && annotation.protospacer_state === "ineligible") {
        wrap.append(node("span", `Protospacer search skipped: ${annotation.skip_reason.replaceAll("_", " ")}`, "match-status"));
      } else if (annotation[kind + "_state"] === "complete") {
        wrap.append(node("span", kind === "database" ? "No exact catalogue match" : "No accepted protospacer match", "match-status"));
      }
      continue;
    }
    const details = node("details", undefined, `match-details ${kind}-evidence`);
    details.append(node("summary", `${count.toLocaleString()} ${label} ${count === 1 ? "match" : "matches"}`));
    const list = node("ul");
    const append = hits => {
      for (const hit of hits) {
        const item = node("li"), url = evidenceUrl(hit.url);
        const link = node(url ? "a" : "span", evidenceLabel(hit, kind));
        if (url) { link.href = url; link.target = "_blank"; link.rel = "noopener noreferrer"; }
        item.append(link);
        if (kind === "protospacer") {
          item.append(node("span", hit.subject_title, "match-description"));
          item.append(node("span", `E-value ${hit.evalue} · bit score ${hit.bitscore}`, "match-description"));
        } else if (hit.certainty_score != null) {
          item.append(node("span", `Source confidence ${Number(hit.certainty_score).toFixed(3)}`, "match-description"));
        }
        list.append(item);
      }
    };
    const initial = annotation[kind + "_matches"] || [];
    append(initial); details.append(list);
    let offset = initial.length;
    if (offset < count && loadMore) {
      const button = node("button", "Show more matches", "secondary more-matches");
      button.type = "button";
      button.addEventListener("click", async () => {
        button.disabled = true;
        try {
          const data = await loadMore(annotation.query_id, kind, offset);
          append(data.matches); offset += data.matches.length;
          if (offset >= data.total) button.remove();
        } catch (error) { onError(error); }
        finally { button.disabled = false; }
      });
      details.append(button);
    } else if (offset < count) {
      details.append(node("p", "Additional matches are available in the complete result download.", "match-status"));
    }
    wrap.append(details);
  }
  return wrap;
}
