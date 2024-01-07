function [] = tima_mapping_recombine_rows(Line_DIR,Nrows,Ncols,GeoRefFileName,newfilename_prefix,top_layer)
%% TIMA_MAPPING_RECOMBINE_ROWS (TI Mars Analogs)
%   Function to generate thermal inertia maps from line-by-line thermal conductivity data.
%
% Description
%   This script is meant to be used to recombine raster data which was output line by line
%   in a parallel task. This function deletes line files.
%
% Input Parameters
%   Line_DIR: Directory of files - TK_Line_%u.txt, Depth_Line_%u.txt, fval_Line_%u.txt (char)
%   Nrows: Number of rows (vector)
%   Ncols: number of columns (vector)
%   GeoRefFileName: Full path and name of Nrows x Ncols geotiff in proper reference (char)
%   newfilename_prefix: Full path and name prefix for geotiffs to be generated (char)
%   top_layer: Thickness in meters of topmost subsurface grid layer used to in thermal inertia derivation (vector)
%
% Outputs:
%   Produces gridded .mat and geotiff raster data of thermal inertia, depth-to-lower-layer, and fval goodness of fit.
%
% Author
%    Ari Koeppel -- Copyright 2023
%
% See also 
%   TIMA_HEAT_TRANSFER TIMA_INITIALIZE TIMA_LATENT_HEAT_MODEL TIMA_LN_PRIOR TIMA_SENSIBLE_HEAT_MODEL TIMA_GWMCMC TIMA_COMBINE_ROWS

% Line_DIR = '/scratch/ak2223/TI_2D_Mapping/Parallel_EE2022/Output_10cm_15Sec/';
% Nrows = 1844; Ncols = 885;
% GeoRefFileName = 'X:\akoeppel\TI_Modeling\Parallel\IcelandEmergingEsker2022_10cm\EmergingEsker_Clipped_10cm_GeoReference.tif';
% newfilename_prefix = 'X:\akoeppel\TI_Modeling\Parallel\IcelandEmergingEsker2022_10cm\Output\EmergingEsker_10cm_15sec_';

LineTKData = NaN([Nrows Ncols]);
LineDepthData = NaN([Nrows Ncols]);
LinefvalData = NaN([Nrows Ncols]);
for rows = 1:Nrows
    LineTKFile = [Line_DIR,sprintf('TK_Line_%u.txt',rows)];
    LineDepthFile = [Line_DIR,sprintf('Depth_Line_%u.txt',rows)];
    LinefvalFile = [Line_DIR,sprintf('fval_Line_%u.txt',rows)];
    if isfile(LineTKFile)
        LineTKData(rows,:) = readmatrix(LineTKFile);
	    LineDepthData(rows,:) = readmatrix(LineDepthFile);
        LinefvalData(rows,:) = readmatrix(LinefvalFile);
    end
end

c = fix(clock);
save([newfilename_prefix,sprintf('%02.0f%02.0f-%s_RAW.mat',c(4),c(5),date)],'LineTKData','LineDepthData','LinefvalData')
LineTKData(LineTKData==0)=NaN;
LineTKData(LineTKData == 0.65535) = NaN;
LineDepthData(LineDepthData<0)=0;
LinefvalData(LinefvalData==0)=NaN;
LinefvalData(isnan(LineTKData)) = NaN;
save([newfilename_prefix,sprintf('%02.0f%02.0f-%s_REFINED.mat',c(4),c(5),date)],'LineTKData','LineDepthData','LinefvalData')

delete([Line_DIR,'*.txt'])

% load('C:\Users\ahkoe\Desktop\NAU_PhD_Work\TI_Earth_2D\EmergingEsker2022\RESULTS_Refined_10cm_15Sec_0858-07-Jan-2024.mat');

%% Place Raster into geotiff

%The top layer is currently always sediment (not ice)
%It should be possible to change the properties of the top layer to water
%and adjust latent heat to reflect standing water.
LineDepthData(LineDepthData<=top_layer) = 0;
LineCPData = ones(size(LinefvalData));
LinerhoData = ones(size(LinefvalData));
LinerhoData(LineDepthData>top_layer) = 1864; %Density measured for sediment GM_01A:1235.13,GM_01B:1373.86,GM_02:1204.40,GM_03:1190.79,GM_04:965.73,GM_05:1223.04,GM_06:914.00,GM_07:1501.17,GM_08:1864.38,GM_09:1896.20,GM_10:1936.59,GM_11:1638.10,GM_12:1430.96
LinerhoData(LineDepthData<=top_layer) = 950; %Density between water and Ice
LineCPData(LineDepthData<=top_layer) = 2100; % Cp Ice (maybe larger with water (4000) and rock (600) mixed in)
LineCPData(LineDepthData>top_layer) = 700; %~Cp sandy(800) basalt(600)
% LineTKData(LineDepthData<=Layer_size_B(1)) = 2; %k ice (maybe smaller with water (0.690) vs pure ice (2.22))

status = copyfile(GeoRefFileName,[newfilename_prefix,'TI.tif'],'f');
t = Tiff(newfilename_prefix,'r+');
T_File = sqrt(LineTKData.*LineCPData.*LinerhoData);
setTag(t,'BitsPerSample',32); %An internal note to program to treat data as 32 Bit float rather than uint16 -not actually written to file
setTag(t,'SampleFormat',Tiff.SampleFormat.IEEEFP); %An internal note to program to treat data as ~float rather than uint16 -not actually written to file
write(t,single(T_File));
close(t);

status = copyfile(GeoRefFileName,[newfilename_prefix,'Depth.tif'],'f');
tt = Tiff(newfilename_prefix,'r+');
T_File = LineDepthData;
setTag(tt,'BitsPerSample',32); %An internal note to program to treat data as 32 Bit float rather than uint16 -not actually written to file
setTag(tt,'SampleFormat',Tiff.SampleFormat.IEEEFP); %An internal note to program to treat data as ~float rather than uint16 -not actually written to file
write(tt,single(T_File));
close(tt);

status = copyfile(GeoRefFileName,[newfilename_prefix,'fval.tif'],'f');
ttt = Tiff(newfilename_prefix,'r+');
T_File = LinefvalData;
setTag(ttt,'BitsPerSample',32); %An internal note to program to treat data as 32 Bit float rather than uint16 -not actually written to file
setTag(ttt,'SampleFormat',Tiff.SampleFormat.IEEEFP); %An internal note to program to treat data as ~float rather than uint16 -not actually written to file
write(ttt,single(T_File));
close(ttt);
end