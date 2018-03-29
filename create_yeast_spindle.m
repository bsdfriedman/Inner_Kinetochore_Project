function create_yeast_spindle

% on/off codes
create_horns = 1; % 0 does not create horns, 1 does
create_ndc80_ring = 1; %0 does not create the ndc80 ring, 1 does

% parameters
num_beads_mt = 8; % number of beads in the MT ring
mt_base_diam_nm = 25; % diameter of the microtubule in nm
mt_horn_diam_nm = 60; % diameter of the horn portion of the microtubule in nm
spindle_outer_diam_nm = 250; % diameter of the outer MTs in the spindle
spindle_inner_diam_nm = 125; % diameter of the outer MTs in the spindle
spindle_length_nm = 1500; % length of the spindle in nm
pericen_length_nm = 800; % length of the pericentromere in nm
inner_mt_overlap_nm = 100; % overlap of the inner sets of MTs in nm