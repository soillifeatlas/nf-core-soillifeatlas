/**
 * command_generator.js
 *
 * Pure function: form state (matching `form_schema.json`) -> `nextflow run`
 * CLI string. No runtime dependencies; works as an ES module and via CommonJS.
 *
 * Parameter names match the real params declared in nextflow.config /
 * nextflow_schema.json of nf-core-soillifeatlas v0.1.0-dev. Defaults mirror
 * those same sources so that a form state at defaults produces the shortest
 * possible command.
 *
 * @typedef {Object} FormState
 * @property {Object} upload
 * @property {string} upload.soil_intensity    - Path to per-sample intensity matrix (parquet/csv).
 * @property {string} upload.soil_mgf          - Path to MS2 spectra MGF.
 * @property {string} [upload.sample_metadata] - Optional TSV/CSV with treatment metadata.
 * @property {Object} instrument
 * @property {'positive'} instrument.ion_mode  - v0.1: only 'positive'.
 * @property {string} instrument.is_compound   - EquiSPLASH IS compound name (e.g. 'LPE_18d7').
 * @property {number} instrument.is_spiked_pmol- IS spike amount in pmol.
 * @property {string} instrument.atlas_version - Atlas release tag.
 * @property {Object} [advanced]
 * @property {number} [advanced.rie_floor]
 * @property {boolean} [advanced.restrict_archaea_to_archlips]
 * @property {string[]} [advanced.decomposition_methods]
 * @property {string} [advanced.primary_method]
 * @property {boolean} [advanced.skip_gnps]
 * @property {boolean} [advanced.skip_annotate]
 */

const DEFAULTS = Object.freeze({
  rie_floor: 0.2,
  restrict_archaea_to_archlips: true,
  decomposition_methods: ['nnls', 'std_bc', 'enriched_bc', 'fc_weighted_bc'],
  primary_method: 'fc_weighted_bc',
  skip_gnps: false,
  skip_annotate: true,
  IS_reference_compound: 'LPE_18d7',
  IS_spiked_pmol: 100,
  atlas_version: 'v0.1.0',
});

const PIPELINE_REVISION = 'v0.1.0-dev';
const PIPELINE_NAME = 'soillifeatlas/nf-core-soillifeatlas';
const DEFAULT_PROFILE = 'local';

/** Shell-quote a value for safe inclusion in the generated command. */
function shq(v) {
  if (v === null || v === undefined) return "''";
  const s = String(v);
  // Avoid quoting simple identifier-like tokens; quote everything else.
  return /^[A-Za-z0-9_./\-+=:,]+$/.test(s) ? s : `'${s.replace(/'/g, "'\\''")}'`;
}

/** Compare two primitive arrays (order-sensitive). */
function arraysEqual(a, b) {
  return Array.isArray(a) && Array.isArray(b) && a.length === b.length &&
    a.every((v, i) => v === b[i]);
}

/**
 * Generate a nextflow run command from form state.
 *
 * Only non-default advanced flags are emitted — a form at defaults produces
 * the minimal "upload + instrument" command.
 *
 * @param {FormState} formState
 * @returns {string} Multi-line `nextflow run ...` command.
 */
export function generateNextflowCommand(formState) {
  const { upload = {}, instrument = {}, advanced = {} } = formState || {};

  const flags = [];
  const push = (flag, val) => flags.push(`  --${flag} ${shq(val)}`);

  // --- required: inputs -----------------------------------------------------
  if (upload.soil_intensity)  push('soil_intensity',  upload.soil_intensity);
  if (upload.soil_mgf)        push('soil_mgf',        upload.soil_mgf);
  if (upload.sample_metadata) push('sample_metadata', upload.sample_metadata);

  // --- instrument / atlas (emit when non-default so reviewers see the intent)
  if (instrument.atlas_version && instrument.atlas_version !== DEFAULTS.atlas_version) {
    push('atlas_version', instrument.atlas_version);
  }
  if (instrument.is_compound && instrument.is_compound !== DEFAULTS.IS_reference_compound) {
    push('IS_reference_compound', instrument.is_compound);
  }
  if (instrument.is_spiked_pmol !== undefined &&
      Number(instrument.is_spiked_pmol) !== DEFAULTS.IS_spiked_pmol) {
    push('IS_spiked_pmol', instrument.is_spiked_pmol);
  }
  // ion_mode = 'positive' maps to enable_neg=false (default); nothing to emit.
  // If we ever support 'negative', emit --enable_neg true here.

  // --- advanced (only diffs from defaults) ----------------------------------
  if (advanced.rie_floor !== undefined &&
      Number(advanced.rie_floor) !== DEFAULTS.rie_floor) {
    push('RIE_floor', advanced.rie_floor);
  }
  if (advanced.restrict_archaea_to_archlips !== undefined &&
      advanced.restrict_archaea_to_archlips !== DEFAULTS.restrict_archaea_to_archlips) {
    push('restrict_archaea_to_archlips', advanced.restrict_archaea_to_archlips);
  }
  if (Array.isArray(advanced.decomposition_methods) &&
      !arraysEqual(advanced.decomposition_methods, DEFAULTS.decomposition_methods)) {
    push('decomposition_methods', advanced.decomposition_methods.join(','));
  }
  if (advanced.primary_method && advanced.primary_method !== DEFAULTS.primary_method) {
    push('primary_method', advanced.primary_method);
  }
  if (advanced.skip_gnps !== undefined && advanced.skip_gnps !== DEFAULTS.skip_gnps) {
    push('skip_gnps', advanced.skip_gnps);
  }
  if (advanced.skip_annotate !== undefined && advanced.skip_annotate !== DEFAULTS.skip_annotate) {
    push('skip_annotate', advanced.skip_annotate);
  }

  const header = [
    `nextflow run ${PIPELINE_NAME} \\`,
    `  -r ${PIPELINE_REVISION} \\`,
    `  -profile ${DEFAULT_PROFILE} \\`,
    `  --mode apply \\`,
    `  --outdir ./results`,
  ];
  return flags.length ? header.join('\n') + ' \\\n' + flags.join(' \\\n') : header.join('\n');
}

export default generateNextflowCommand;

// CommonJS interop (Node + bundlers both).
if (typeof module !== 'undefined' && module.exports) {
  module.exports = generateNextflowCommand;
  module.exports.default = generateNextflowCommand;
  module.exports.generateNextflowCommand = generateNextflowCommand;
}

/* example

// --- Test case 1: minimal (defaults everywhere) ----------------------------
const demo1 = {
  upload: { soil_intensity: 'soil.parquet', soil_mgf: 'soil.mgf' },
  instrument: { ion_mode: 'positive', is_compound: 'LPE_18d7',
                is_spiked_pmol: 100, atlas_version: 'v0.1.0' },
};
// generateNextflowCommand(demo1) ->
// nextflow run soillifeatlas/nf-core-soillifeatlas \
//   -r v0.1.0-dev \
//   -profile local \
//   --mode apply \
//   --outdir ./results \
//   --soil_intensity soil.parquet \
//   --soil_mgf soil.mgf

// --- Test case 2: with metadata + a tweaked RIE floor ----------------------
const demo2 = {
  upload: {
    soil_intensity: '/data/climgrass/soil.parquet',
    soil_mgf: '/data/climgrass/soil.mgf',
    sample_metadata: '/data/climgrass/meta.tsv',
  },
  instrument: { ion_mode: 'positive', is_compound: 'LPE_18d7',
                is_spiked_pmol: 100, atlas_version: 'v0.1.0' },
  advanced: { rie_floor: 0.15 },
};
// generateNextflowCommand(demo2) ->
// nextflow run soillifeatlas/nf-core-soillifeatlas \
//   -r v0.1.0-dev \
//   -profile local \
//   --mode apply \
//   --outdir ./results \
//   --soil_intensity /data/climgrass/soil.parquet \
//   --soil_mgf /data/climgrass/soil.mgf \
//   --sample_metadata /data/climgrass/meta.tsv \
//   --RIE_floor 0.15

// --- Test case 3: fast mode (2 methods, skip GNPS, switch primary) ---------
const demo3 = {
  upload: { soil_intensity: 'a.parquet', soil_mgf: 'a.mgf' },
  instrument: { ion_mode: 'positive', is_compound: 'PC_15-18d7',
                is_spiked_pmol: 50, atlas_version: 'v0.1.0' },
  advanced: {
    decomposition_methods: ['nnls', 'fc_weighted_bc'],
    primary_method: 'nnls',
    skip_gnps: true,
    restrict_archaea_to_archlips: false,
  },
};
// generateNextflowCommand(demo3) ->
// nextflow run soillifeatlas/nf-core-soillifeatlas \
//   -r v0.1.0-dev \
//   -profile local \
//   --mode apply \
//   --outdir ./results \
//   --soil_intensity a.parquet \
//   --soil_mgf a.mgf \
//   --IS_reference_compound PC_15-18d7 \
//   --IS_spiked_pmol 50 \
//   --restrict_archaea_to_archlips false \
//   --decomposition_methods nnls,fc_weighted_bc \
//   --primary_method nnls \
//   --skip_gnps true

*/
