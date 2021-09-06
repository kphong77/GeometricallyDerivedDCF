
close all; clear all; clc

% the number of sampled points in a radial projection
Nsamples = 64       
% double sampling in MRI (actual sampling in a radial projection)
NM = Nsamples*2     
% the number of radial projections
Nproj = 100            
% the number of repeats
Nrepeat = 1         
% select KWIC function
sel_WC = 0          
% projection angle type for 2D
%    0=regular angle (evenly-spaced angle)
%    otherwise=tiny golden angle 
sel_Angle2D = 0   


clear  AngleInfo
switch sel_Angle2D
    case 0   % regular angle
        RA = pi/(Nproj*Nrepeat)
        AngleInfo = ((0:Nproj*Nrepeat-1)*RA)';
    otherwise % tiny golden angle
        Ntiny = sel_Angle2D;   
        TGA = pi/((sqrt(5)+1)/2 + (Ntiny-1) )
        AngleInfo = mod( (0:Nproj*Nrepeat-1)*TGA, 2*pi )';
end
AngleInfo = reshape( AngleInfo, [Nproj Nrepeat] );
size(AngleInfo)

Option.Nsamples = NM;
Option.AnglePrecision = 0;   % full precision
Option.Display = 1;
if sel_WC == 0
    Option.WeightedContrast = 0
else
    LocOfMainContrastRay = [3 5 7];   % selected projections in this DEMO setting
    Option.WeightedContrast = ones( NM, Nproj, Nrepeat );
    Option.WeightedContrast( NM/2+1-2:NM/2+1+2, :, : ) = 0;
    for ii = 1:Nrepeat
        Option.WeightedContrast(  NM/2+1,              LocOfMainContrastRay(ii),                              ii ) = 1;
        Option.WeightedContrast( [NM/2+1-1  NM/2+1+1], LocOfMainContrastRay(ii)-1:LocOfMainContrastRay(ii)+1, ii ) = 1;
        Option.WeightedContrast( [NM/2+1-2  NM/2+1+2], LocOfMainContrastRay(ii)-2:LocOfMainContrastRay(ii)+2, ii ) = 1;
    end
end

DCF = OptimalDCF_2DRadial( AngleInfo, Option );

figure(1); imagesc( DCF(:,:,1), [0 1] ); colormap jet
axis image
yticks( [1  NM/2+1  NM] )
yticklabels( { num2str(NM/2), '0', num2str(NM/2-1) } )
set(gca, 'fontsize', 30, 'fontweight', 'bold' )
ylabel( 'Radius' )
xlabel( 'Projection' )
colorbar

