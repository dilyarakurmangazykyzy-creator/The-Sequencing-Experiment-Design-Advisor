"use strict";

const $ = id => document.getElementById(id);
const form = $("design-form");
const percentFields = new Set(["allele_fraction", "detection_probability", "mapping_fraction", "usable_fraction"]);
const stringFields = new Set(["organism", "application", "platform", "procurement", "assembly_goal", "expression_goal"]);
const money = value => new Intl.NumberFormat("en-US", {style:"currency", currency:"USD", maximumFractionDigits:0}).format(value);
const moneyExact = value => new Intl.NumberFormat("en-US", {style:"currency", currency:"USD", minimumFractionDigits:2, maximumFractionDigits:2}).format(value);
const number = (value, digits = 1) => new Intl.NumberFormat("en-US", {maximumFractionDigits:digits}).format(value);
const percentage = (value, digits = 1) => number(value * 100, digits) + "%";
const escapeHTML = value => String(value ?? "").replace(/[&<>"']/g, char => ({"&":"&amp;","<":"&lt;",">":"&gt;",'"':"&quot;","'":"&#39;"}[char]));
const finite = value => typeof value === "number" && Number.isFinite(value);
let configBundle = null;
let currentResult = null;
let requestNumber = 0;
let toastTimer = null;
const costColors = ["#123f37", "#487f63", "#8eaf75", "#c4dfa2", "#287c91", "#79a9b6", "#d8c59a"];
const costLabels = {library_preparation:"Library preparation", sample_qc:"Sample QC", sequencing_consumables:"Sequencing consumables", instrument_time:"Instrument allocation", cloud_compute:"Analysis compute", storage:"Storage", contingency:"Contingency"};
const platformColors = {illumina:"#287c91", nanopore:"#bf6b39", pacbio:"#749842"};

async function fetchJSON(url, options = {}) {
  const response = await fetch(url, options);
  const value = await response.json();
  if (!response.ok) throw new Error(value.error || `Request failed (${response.status}).`);
  return value;
}

function toast(message) {
  $("toast").textContent = message;
  $("toast").hidden = false;
  clearTimeout(toastTimer);
  toastTimer = setTimeout(() => $("toast").hidden = true, 3500);
}

function downloadJSON(value, name) {
  const blob = new Blob([JSON.stringify(value, null, 2) + "\n"], {type:"application/json"});
  const url = URL.createObjectURL(blob);
  const link = document.createElement("a");
  link.href = url;
  link.download = name;
  link.click();
  setTimeout(() => URL.revokeObjectURL(url), 1000);
  toast("JSON saved. The assumptions remain attached to the design.");
}

function applyScenario(input) {
  const values = {...configBundle.defaults, ...input};
  if (input.organism && input.genome_size_mb === undefined) {
    values.genome_size_mb = {bacteria:4.6, yeast:12.1, human:3100, custom:100}[input.organism];
  }
  for (const [name, value] of Object.entries(values)) {
    const control = $(name);
    if (control) control.value = percentFields.has(name) ? Number((value * 100).toFixed(8)) : value;
  }
  updateQuestionFields();
}

function readScenario() {
  const input = {};
  for (const [key, value] of new FormData(form).entries()) {
    if (stringFields.has(key)) input[key] = value;
    else input[key] = Number(value) / (percentFields.has(key) ? 100 : 1);
  }
  return input;
}

function updateQuestionFields() {
  const application = $("application").value;
  for (const name of ["variant", "expression", "assembly"]) $(name + "-fields").hidden = application !== name;
  const isoforms = $("expression_goal").value === "full_length_isoforms";
  $("pairs-field").hidden = isoforms;
  $("long-count-field").hidden = !isoforms;
}

async function runDesign() {
  if (!form.reportValidity()) return;
  const generation = ++requestNumber;
  $("form-error").hidden = true;
  $("run-button").disabled = true;
  try {
    const result = await fetchJSON("/api/recommend", {method:"POST", headers:{"Content-Type":"application/json"}, body:JSON.stringify(readScenario())});
    if (generation !== requestNumber) return;
    currentResult = result;
    renderResult(result);
    $("loading").hidden = true;
    $("result-content").hidden = false;
  } catch (error) {
    if (generation !== requestNumber) return;
    $("form-error").textContent = error.message;
    $("form-error").hidden = false;
    $("loading").textContent = "The design could not be calculated. Review the input message.";
  } finally {
    if (generation === requestNumber) $("run-button").disabled = false;
  }
}

function renderResult(result) {
  const r = result.recommendation;
  const input = result.input;
  const cost = r.cost;
  const detection = r.detection;
  const statusLabels = {within_estimate:"FITS ESTIMATE", over_budget:"OVER BUDGET", incomplete_design:"INCOMPLETE DESIGN", unsupported:"UNSUPPORTED OBJECTIVE"};
  $("status-pill").textContent = statusLabels[result.status] || result.status;
  $("status-pill").className = "status-pill" + (result.status === "within_estimate" ? "" : " bad");
  $("design-platform").textContent = r.label;
  $("design-instrument").textContent = r.instrument;
  $("design-read-length").textContent = r.read_length;
  $("cost-value").textContent = money(cost.total_usd);
  $("cost-note").textContent = r.within_budget ? `${money(Math.max(0,input.budget_usd - cost.total_usd))} below budget · assumption-based` : `${money(r.budget_gap_usd)} budget gap`;
  $("depth-label").textContent = input.application === "expression" ? "Assigned fragments / sample" : "Effective mean depth";
  $("depth-value").textContent = input.application === "expression" ? number(r.effective_fragments_per_sample / 1e6, 1) + "M" : number(r.effective_depth_x, 0) + "×";
  $("depth-note").textContent = input.application === "expression" ? `${input.sample_count} biological samples · no genomic depth` : `${number(r.nominal_depth_x, 1)}× planned raw nominal depth`;
  if (detection.type === "idealized_allele_sampling") {
    $("limit-label").textContent = "Ideal sampling threshold";
    $("limit-value").textContent = percentage(detection.minimum_allele_fraction_at_target_probability, 1);
    $("limit-note").textContent = `At ${percentage(detection.requested_probability, 1)} chance of ≥${detection.minimum_alt_reads} alt reads · specificity unvalidated`;
  } else if (detection.type === "idealized_transcript_sampling") {
    $("limit-label").textContent = "Fragment fraction threshold";
    $("limit-value").textContent = number(detection.minimum_fragment_fraction_ppm, 2) + " ppm";
    $("limit-note").textContent = `≥${detection.minimum_count_threshold} counts at ${percentage(detection.requested_probability)} chance · DE power unvalidated`;
  } else {
    $("limit-label").textContent = "Required span + anchors";
    $("limit-value").textContent = number(detection.required_span_bp / 1000, 1) + " kb";
    $("limit-note").textContent = "Length filter · completeness unvalidated";
  }
  $("design-headline").textContent = result.headline + ". " + (r.feasible ? "Confirm library-specific performance before ordering." : "Review the gap or objective before using this design.");
  $("reasons").innerHTML = r.reasons.map(text => `<li>${escapeHTML(text)}</li>`).join("");
  $("caveats").innerHTML = r.caveats.length ? r.caveats.map(text => `<p>${escapeHTML(text)}</p>`).join("") : "Planning depth is provisional. The empirical endpoint has not validated variant sensitivity.";
  renderCost(cost);
  $("alternatives-body").innerHTML = result.alternatives.map(candidate => {
    const chosen = candidate.platform === r.platform;
    const decision = !candidate.compatible ? "Objective unsupported" : !candidate.biological_design_adequate ? "Sample design incomplete" : !candidate.within_budget ? "Over budget" : "Fits estimate";
    const required = !candidate.compatible ? "Not modeled for this goal" : candidate.effective_depth_x === null ? number(candidate.effective_fragments_per_sample / 1e6) + "M fragments" : number(candidate.effective_depth_x, 0) + "× effective";
    return `<tr><td>${escapeHTML(candidate.label)}<small>${escapeHTML(candidate.read_length)}</small></td><td>${required}</td><td>${candidate.compatible ? money(candidate.cost.total_usd) : "—"}${candidate.compatible && !candidate.within_budget ? `<small>Gap ${money(candidate.budget_gap_usd)}</small>` : ""}</td><td><span class="choice${!candidate.feasible ? " bad" : chosen ? " chosen" : ""}">${chosen && candidate.feasible ? "Selected" : decision}</span></td></tr>`;
  }).join("");
  renderEvidence(result.evidence);
  renderMethod(result);
  renderSources(result);
}

function renderCost(cost) {
  const entries = Object.entries(cost.components_usd);
  const total = entries.reduce((sum, item) => sum + item[1], 0);
  let offset = 0;
  const rects = entries.map(([key, value], index) => {
    const width = total ? value / total * 1000 : 0;
    const rect = `<rect x="${offset}" y="0" width="${width}" height="22" fill="${costColors[index]}"><title>${costLabels[key]}: ${money(value)}</title></rect>`;
    offset += width;
    return rect;
  });
  $("cost-chart").innerHTML = `<svg viewBox="0 0 1000 22" preserveAspectRatio="none" role="img" aria-label="Project cost decomposition">${rects.join("")}</svg>`;
  $("cost-legend").innerHTML = entries.map(([key,value],index) => `<div><svg class="cost-key" viewBox="0 0 8 8" aria-hidden="true"><rect width="8" height="8" fill="${costColors[index]}"/></svg><span>${costLabels[key]}</span><strong>${moneyExact(value)}</strong></div>`).join("");
  $("cost-range").textContent = `Sensitivity range ${money(cost.planning_range_usd[0])}–${money(cost.planning_range_usd[1])}`;
  $("raw-volume").textContent = `${number(cost.raw_gb, 3)} Gb raw project yield`;
  $("cost-provenance").textContent = cost.provenance;
  $("cost-procurement").textContent = cost.procurement === "dedicated"
    ? `${cost.dedicated_flow_cells} whole flow cell(s) · ${number(cost.billed_gb, 2)} Gb capacity · ${number(cost.allocated_instrument_hours, 1)} instrument hours. Unused capacity is charged. Capacity is a planning estimate, not guaranteed usable yield.`
    : `Shared allocation: ${number(cost.flow_cell_equivalents, 5)} flow-cell equivalent(s), ${number(cost.allocated_instrument_hours, 3)} instrument hours; includes a minimum sequencing allocation. Actual core scheduling and minimum orders vary.`;
}

function chart(series, {title, xLabel="Nominal sampled depth (×)", target=null, maxX=null} = {}) {
  const all = series.flatMap(item => item.points).filter(point => finite(point.x) && finite(point.y));
  if (!all.length) return "<p class='help-note'>No measured points are available for this endpoint.</p>";
  const width=550,height=255,left=48,right=15,top=21,bottom=44;
  const xMaximum = maxX || Math.max(1, ...all.map(point => point.x));
  const x = value => left + value / xMaximum * (width-left-right);
  const y = value => height-bottom - Math.max(0,Math.min(1,value)) * (height-top-bottom);
  const nodes = [];
  for (let tick=0;tick<=4;tick++) {
    const value=tick/4;
    nodes.push(`<line x1="${left}" y1="${y(value)}" x2="${width-right}" y2="${y(value)}" class="chart-gridline"/><text x="${left-9}" y="${y(value)+3}" text-anchor="end" class="chart-axis">${number(value*100,0)}%</text>`);
  }
  for (let tick=0;tick<=4;tick++) {
    const value=tick/4*xMaximum;
    nodes.push(`<text x="${x(value)}" y="${height-bottom+20}" text-anchor="middle" class="chart-axis">${number(value,xMaximum<5?2:0)}</text>`);
  }
  if (finite(target)) nodes.push(`<line x1="${left}" y1="${y(target)}" x2="${width-right}" y2="${y(target)}" class="chart-target"/><text x="${width-right}" y="${y(target)-6}" text-anchor="end" class="chart-title">target ${percentage(target)}</text>`);
  for (const item of series) {
    const points = item.points.filter(point=>finite(point.x)&&finite(point.y)).sort((a,b)=>a.x-b.x);
    if (!points.length) continue;
    if (points.every(point=>finite(point.low)&&finite(point.high))) {
      const bounds = points.map(point=>`${x(point.x)},${y(point.low)}`).concat([...points].reverse().map(point=>`${x(point.x)},${y(point.high)}`)).join(" ");
      nodes.push(`<polygon points="${bounds}" fill="${item.color}" opacity="0.13"/>`);
    }
    nodes.push(`<polyline points="${points.map(point=>`${x(point.x)},${y(point.y)}`).join(" ")}" stroke="${item.color}" stroke-width="2.6" fill="none"/>`);
    for (const point of points) nodes.push(`<circle cx="${x(point.x)}" cy="${y(point.y)}" r="3.2" fill="${item.color}"><title>${escapeHTML(item.label)} · ${number(point.x,3)}× · ${percentage(point.y,2)}</title></circle>`);
  }
  nodes.push(`<text x="${(left+width-right)/2}" y="${height-7}" text-anchor="middle" class="chart-title">${escapeHTML(xLabel)}</text>`);
  return `<svg xmlns="http://www.w3.org/2000/svg" viewBox="0 0 ${width} ${height}" role="img" aria-label="${escapeHTML(title)}">${nodes.join("")}</svg>`;
}

function renderEvidence(evidence) {
  const saturation = evidence.saturation;
  const supported = ["anchor_supported_reference_breadth", "reference_breadth_primary_alignments_mapq_ge_20"];
  const records = saturation && supported.includes(saturation.endpoint) && typeof saturation.platforms === "object" ? saturation.platforms : {};
  const series = Object.entries(records || {}).filter(([,dataset])=>dataset&&Array.isArray(dataset.points)).map(([platform,dataset]) => ({platform, dataset:{...dataset,points:dataset.points.map(point=>({...point,nominal_depth_x:point.nominal_depth_x_mean ?? point.nominal_depth_x}))}, label:platform === "illumina" ? "Illumina" : platform === "nanopore" ? "Nanopore" : platform, color:platformColors[platform]||"#749842"}));
  const hasPoints = series.some(item=>item.dataset.points.some(point=>finite(point.nominal_depth_x)&&finite(point.coverage_1x_mean)));
  $("evidence-status").textContent = hasPoints ? "LIMITED BACTERIAL PILOT" : "EMPIRICAL EVIDENCE INCOMPLETE";
  const standard = saturation?.endpoint === "reference_breadth_primary_alignments_mapq_ge_20";
  $("evidence-scope").textContent = hasPoints ? (standard ? "Primary alignments with MAPQ ≥20 support the measured reference breadth. Only this bacterial mapping endpoint was measured; biological conclusions require their own validation." : "A custom anchor mapper estimates support for reference positions in one bacterial pilot. This is algorithm-specific evidence; it is not a standard aligner benchmark.") : evidence.issues.join(" ");
  const mappingCaption = $("mapping-chart").parentElement.querySelector("figcaption p");
  $("mapping-chart").parentElement.querySelector("figcaption h3").textContent = standard ? "Primary read mapping rate" : "Read anchoring rate";
  mappingCaption.textContent = standard ? "Fraction of sampled reads mapped by the pilot aligner" : "Algorithm-specific; not a standard aligner mapping rate";
  const breadthCaption = $("breadth-chart").parentElement.querySelector("figcaption p");
  breadthCaption.textContent = standard ? "Reference fraction with ≥1 primary aligned read (MAPQ ≥20)" : "Fraction supported by at least one anchored read";
  $("evidence-empty").hidden = hasPoints;
  $("evidence-charts").hidden = !hasPoints;
  const maximumDepth = Math.max(1,...series.flatMap(item=>item.dataset.points).map(point=>point.nominal_depth_x).filter(finite));
  if (hasPoints) {
    $("breadth-chart").innerHTML = chart(series.map(item=>({...item,points:item.dataset.points.map(point=>({x:point.nominal_depth_x,y:point.coverage_1x_mean,low:point.coverage_1x_min,high:point.coverage_1x_max}))})),{title:"Measured bacterial reference breadth against nominal sampled depth",maxX:maximumDepth});
    $("mapping-chart").innerHTML = chart(series.map(item=>({...item,points:item.dataset.points.map(point=>({x:point.nominal_depth_x,y:point.mapping_rate_mean,low:point.mapping_rate_min,high:point.mapping_rate_max}))})),{title:standard ? "Primary read mapping rate against nominal sampled depth" : "Read anchoring rate against nominal sampled depth",maxX:maximumDepth});
  }
  $("pilot-findings").innerHTML = series.map(item=>{
    const plateau = item.dataset.plateau || {};
    const valid = item.dataset.points.filter(point=>finite(point.nominal_depth_x)&&finite(point.coverage_1x_mean));
    const final = valid.length ? valid.reduce((a,b)=>a.nominal_depth_x>b.nominal_depth_x?a:b) : null;
    let finding = plateau.detected && finite(plateau.depth_x) ? `Candidate breadth plateau at ${number(plateau.depth_x,2)}× nominal sampled depth.` : "No plateau meeting the declared rule was established in the sampled range.";
    if (final) finding += ` Highest sampled depth: ${number(final.nominal_depth_x,2)}×; ≥1× breadth: ${percentage(final.coverage_1x_mean,2)}.`;
    const method = item.dataset.method || item.dataset.endpoint || saturation?.endpoint || "";
    const metadata = item.dataset.metadata || {};
    const accession = metadata.run_accession || metadata.accession || "Accession recorded in the run manifest";
    const organism = metadata.scientific_name || "Bacterial pilot";
    const criterion = (plateau.criterion || "").replace(/At least95%/g, "At least 95%").replace(/breadth>=1x/g, "breadth at ≥1×").replace(/gains<1/g, "gains below 1");
    return `<div class="pilot-finding"><b>${escapeHTML(item.label)} · observed endpoint</b><p>${escapeHTML(accession)} · ${escapeHTML(organism)}</p><p>${escapeHTML(finding)}</p><p>${escapeHTML(criterion)}</p><details><summary>Endpoint &amp; alignment policy</summary><p>${escapeHTML(typeof method === 'string' ? method : JSON.stringify(method))}</p><p>${escapeHTML(item.dataset.endpoint || "")}</p><p>${escapeHTML(item.dataset.caveat || "")}</p><p>${escapeHTML(item.dataset.mapping_rate_depth_caveat || "")}</p></details></div>`;
  }).join("");
}

function renderMethod(result) {
  const input = result.input;
  const r = result.recommendation;
  $("sampling-explanation").textContent = input.application === "variant" ? `At a target allele fraction of ${percentage(input.allele_fraction)}, the model finds the depth needed for at least ${input.min_alt_reads} alternate observations with ${percentage(input.detection_probability)} probability. This is a theoretical sampling calculation.` : input.application === "expression" ? "RNA sampling is modeled using assigned fragments and a minimum count threshold. The displayed allele equation explains the general Poisson approach; differential-expression power needs a dispersion model and biological replicates." : "The assembly design applies a provisional effective-depth floor and a repeat-span filter. The allele equation below is not an assembly-quality model.";
  $("sampling-chart").hidden = input.application !== "variant";
  if (input.application === "variant") {
    const maximum = Math.max(60,(r.effective_depth_x || 30)*1.5);
    const points = result.detection_curve.filter(point=>point.effective_depth_x<=maximum).map(point=>({x:point.effective_depth_x,y:point.sampling_probability}));
    $("sampling-chart").innerHTML = chart([{label:"Ideal allele sampling",color:"#4f8a5b",points}],{title:"Theoretical alternate-read sampling probability, not measured recall",xLabel:"Effective mean depth (×)",target:input.detection_probability,maxX:Math.max(...points.map(point=>point.x))});
  }
  $("sampling-details").textContent = input.application === "variant" ? `Selected depth gives ${percentage(r.detection.predicted_sampling_probability,2)} ideal sampling probability. Budget alone supports up to ${number(r.affordable_effective_depth_x,1)}× effective depth under the same cost assumptions. This ceiling is not a recommendation to under-sequence.` : r.detection.interpretation;
  $("assumption-values").textContent = `Assumed mapping ${percentage(input.mapping_fraction)} × post-mapping retention ${percentage(input.usable_fraction)} = ${percentage(result.methodology.effective_fraction)} effective fraction. Genome: ${number(input.genome_size_mb,2)} Mb; samples: ${input.sample_count}.`;
  const config = configBundle.config;
  $("assumption-table").innerHTML = `<table><thead><tr><th>Platform</th><th>DNA output / cell</th><th>Shared sequencing</th><th>DNA prep / sample</th></tr></thead><tbody>${Object.values(config.platforms).map(p=>`<tr><td>${escapeHTML(p.label)}</td><td>${number(p.run_yield_gb,0)} Gb</td><td>$${number(p.shared_sequencing_usd_per_gb,2)} / Gb</td><td>${money(p.library_per_sample_usd)}</td></tr>`).join("")}</tbody></table>`;
  $("model-limits").innerHTML = result.limitations.map(text=>`<li>${escapeHTML(text)}</li>`).join("");
}

function renderSources(result) {
  $("benchmark-checks").innerHTML = result.cost_validation.map(item=>{
    let title, detail;
    if (item.kind === "independent_sequencing_service") {
      title = item.name;
      const rates = item.service_usd_per_gb_range;
      const differences = item.relative_difference_pct_range;
      detail = `Published service $${number(rates[0],3)}–$${number(rates[1],3)}/Gb; model sequencing + instrument $${number(item.model_sequencing_plus_instrument_usd_per_gb,3)}/Gb. Relative difference ${number(differences[0],1)}% to ${number(differences[1],1)}% · ${item.comparison_status.replaceAll("_"," ")}.`;
    } else if (item.kind === "independent_preparation_service") {
      title = item.name;
      detail = `Published ${money(item.price_usd)}/sample; model ${money(item.model_per_sample_usd)}/sample. Relative difference ${number(item.relative_difference_pct,1)}%. Project totals may substantially underestimate a service order.`;
    } else {
      title = `${item.platform} reagent arithmetic: ${money(item.published_reagent_usd)} / ${number(item.raw_gb,0)} Gb = $${number(item.derived_usd_per_gb,2)}/Gb`;
      detail = `Model reagent unit difference ${number(item.arithmetic_relative_difference_pct || 0,2)}%. This is an arithmetic check of the model input.`;
    }
    return `<div class="benchmark"><b>${escapeHTML(title)}</b><p>${escapeHTML(detail)}</p><p>${escapeHTML(item.scope)}</p></div>`;
  }).join("");
  $("sources").innerHTML = result.sources.map(source=>{
    let safeURL = "#";
    try { const url = new URL(source.url); if (url.protocol === "https:") safeURL = url.href; } catch (_) {}
    return `<div class="source-link"><a href="${escapeHTML(safeURL)}" target="_blank" rel="noopener noreferrer">${escapeHTML(source.title)} ↗</a><p>${escapeHTML(source.supports)}</p></div>`;
  }).join("");
  $("source-date").textContent = `Source check: ${configBundle.config.checked_on}. Model ${result.model_version} · configuration ${result.config_version}. Export preserves the configuration SHA-256.`;
}

form.addEventListener("submit", event=>{event.preventDefault();runDesign();});
$("organism").addEventListener("change",()=>{
  if ($("organism").value !== "custom") $("genome_size_mb").value = {bacteria:4.6,yeast:12.1,human:3100}[$("organism").value];
});
form.addEventListener("change",()=>{updateQuestionFields();document.querySelectorAll(".preset").forEach(button=>button.classList.remove("active"));runDesign();});
$("export-result").addEventListener("click",()=>currentResult&&downloadJSON(currentResult,"sequencing-design.json"));
$("export-scenario").addEventListener("click",()=>downloadJSON(readScenario(),"sequencing-scenario.json"));
$("export-costs").addEventListener("click",async()=>{
  try { const bundle = await fetchJSON("/api/config");downloadJSON(bundle.config,"sequencing-cost-assumptions.json"); } catch(error) {toast(error.message);}
});
document.querySelectorAll("[role=tab]").forEach((button,index,buttons)=>{
  button.addEventListener("click",()=>{
    buttons.forEach(item=>{const selected=item===button;item.setAttribute("aria-selected",String(selected));item.tabIndex=selected?0:-1;$(item.getAttribute("aria-controls")).hidden=!selected;});
  });
  button.addEventListener("keydown",event=>{
    if (!["ArrowRight","ArrowLeft","Home","End"].includes(event.key)) return;
    event.preventDefault();
    const next=event.key==="Home"?0:event.key==="End"?buttons.length-1:(index+(event.key==="ArrowRight"?1:-1)+buttons.length)%buttons.length;
    buttons[next].click();buttons[next].focus();
  });
  button.tabIndex=index===0?0:-1;
});

(async function initialize() {
  try {
    configBundle = await fetchJSON("/api/config");
    applyScenario({});
    for (const scenario of configBundle.scenarios) {
      const button=document.createElement("button");button.type="button";button.className="preset";button.textContent=scenario.label;
      button.addEventListener("click",()=>{document.querySelectorAll(".preset").forEach(item=>item.classList.remove("active"));button.classList.add("active");applyScenario(scenario.input);runDesign();});
      $("presets").appendChild(button);
    }
    await runDesign();
  } catch(error) {
    $("loading").textContent="Could not load the advisor: "+error.message+". Start the bundled local application and reload.";
    $("run-button").disabled=true;
  }
})();
