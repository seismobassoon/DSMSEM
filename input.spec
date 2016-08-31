# -*- mode: perl -*-
run_name = "Test_name";

# duration of the run
sim_time = 0.5;
mesh_file = "mesh4spec"; # input mesh file
mat_file = "material.input";
dim=3;

snapshots {
    save_snap = false;
    snap_interval = 0.01;
    deselect all;
    select box = -100 -100  100 500 500 150;
    select box = -100 100  -100 500 150 500;
    select box = 100 -100  -100 150 500 500;
};

# Description des capteurs
save_traces = true;
traces_format=hdf5;

capteurs "src" {
    type = single;
    point0 = 25. 25. 25.;
    period = 1;
};

capteurs "x1" {
    type = single;
    point0 = 125. 25. 25.;
    period = 1;
};

capteurs "x2" {
    type = single;
    point0 = 125. 75. 25.;
    period = 1;
};

capteurs "y1" {
    type = single;
    point0 = 25. 225. 25.;
    period = 1;
};

capteurs "y2" {
    type = single;
    point0 = 75. 225. 25.;
    period = 1;
};

capteurs "z1" {
    type = single;
    point0 = 25. 25. 325.;
    period = 1;
};

capteurs "z2" {
    type = single;
    point0 = 75. 75. 325.;
    period = 1;
};

capteurs "xyz" {
    type = single;
    point0 = 198.20508 198.20508 198.20508; # 25 + 300 / sqrt ( 3 )
    period = 1;
};

# Fichier protection reprise
prorep=false;
prorep_iter=1000;
restart_iter=0;


# introduce a source
source {
    # coordinates of the sources ((x,y,z) or (lat,long,R) if rotundity is considered)
    coords = 25.0 25.0 25.0;
    # the numbers before the labels are here to help convert from previous input.spec format
    # Type (1.Impulse, 2.moment Tensor, 3.fluidpulse)
    type = impulse;
    # Direction 0.x,1.y ou 2.z (only for Impulse)
    dir = 1. 0. 0.;
    # Function 1.gaussian,2.ricker,3.tf_heaviside,4.gabor,5.file,6.spice_bench,7.sinus
    func = ricker;
    tau = 0.4;
    freq = 3.;   # source main frequency / cutoff frequency
};

time_scheme {
    accel_scheme = false;  # Acceleration scheme for Newmark
    veloc_scheme = true;   # Velocity scheme for Newmark
    alpha = 0.5;           # alpha (Newmark parameter)
    beta =  0.5;           # beta (Newmark parameter)
    gamma = 1;             # gamma (Newmark parameter)
    courant=0.2;
};

amortissement {
    nsolids = 0;           # number of solids for attenuation (0 if no attenuation)
    atn_band = 10  0.05;   # attenuation period band
    atn_period = 0.2;      # model period 
};
