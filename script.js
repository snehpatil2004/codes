const codes = {

exp7_ofdm: `
PASTE CODE FROM:
ofdm_waveform with diff modulation_exp7
`,

exp8_cu_du: `
PASTE CODE FROM:
cu_du_split_5g_exp8
`,

exp9_numerology: `
PASTE CODE FROM:
exp9_diff5g_numerology
`,

exp10_csi_rs: `
PASTE CODE FROM:
exp10_csi_rs
`,

exp11_prach: `
PASTE CODE FROM:
exp11_prach
`,

exp12_mimo: `
PASTE CODE FROM:
exp12_mimo
`,

evolution_modulation: `
PASTE CODE FROM:
evolution of modulation schemes_6th
`

};

function copyCode(key) {
  navigator.clipboard.writeText(codes[key]);
  alert("Experiment code copied!");
}
