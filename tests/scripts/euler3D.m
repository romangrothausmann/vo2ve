#!/usr/bin/octave -q
##script to execute imEuler3d from http://in.mathworks.com/matlabcentral/fileexchange/33690-geometric-measures-in-2d-3d-images


addpath ("~/octave/imMinkowski/");

arg_list = argv ();
if nargin != 1
  fprintf(stderr, "Usage: %s <input3D.mha>\n", program_name);
  fprintf(stderr, "Decompressing MHAs/MHDs is very slow!\n");
  exit(1)
else
  fprintf(stderr, "Reading 3D data from %s...\n", arg_list{1});
  i3d= mha_read_volume(arg_list{1});#from ~/octave/functions/
  # if nargin == 2
  #   if (arg_list{2} == "-q")
  #     quiet= 1;
  #   endif
  # endif
endif

#printf("dim= %d\n", ndims(i3d))
printf("EPC(6)= %d\n", imEuler3d(i3d, 6))
printf("EPC(26)= %d\n", imEuler3d(i3d, 26))



