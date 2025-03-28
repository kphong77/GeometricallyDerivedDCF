
close all; clear all; clc

% the number of sampled points in a radial projection
Nsamples = 64       
% double sampling in MRI (actual sampling in a radial projection)
NM = Nsamples*2     
% the number of radial projections
if 1
    Nproj = 100            
else
    Nproj = round( (4*pi*(Nsamples)^2)/2 ); % Nyquist sampling 
end
% the number of repeats
Nrepeat = 1        
% KWIC function
sel_WC = 0          
% projection angle type for 3D
%  1=spiral phyllotaxis (Piccini D, et al. MRM 2011;66:1049-1056) 
%  2=golden means (Chan RW, et al. MRM 2009;61:354-363)
sel_Angle3D = 1   


% generate angles
AngleInfo = [];
switch sel_Angle3D
    case 1  % spiral phyllotaxis
        GA = 137.51;  %deg
        AngInc = (2*pi/360)*GA;
        N = Nproj*Nrepeat;
        counter = 0;
        for Interleave = 1:Nrepeat
            for Projection = 1:Nproj
                AngleInfo(Projection,1,Interleave) = mod( counter*AngInc, 2*pi );       %azimuthal angle
                AngleInfo(Projection,2,Interleave) = mod( pi/2*sqrt(counter/N), pi );   %polar angle
                counter = counter + 1;
            end
        end
    case 2  % golden means  
        AngInc_azimuth = 0.4656;
        AngInc_polar = 0.6823;
        AngleInfo(:,1) = mod( (0:Nproj*Nrepeat-1)*AngInc_azimuth, 2*pi );
        AngleInfo(:,2) = mod( (0:Nproj*Nrepeat-1)*AngInc_polar, pi );
        AngleInfo = reshape( AngleInfo, [Nproj 2 Nrepeat] );
end    
size(AngleInfo)


% calculate DCF
Option.Nsamples = NM;
Option.AnglePrecision = 0;   % full precision
Option.WeightedContrast = 0;
Option.Display = 1;

gDCF = gDCF_3DRadial( AngleInfo, Option );
size(gDCF)

figure(1); imagesc( gDCF(:,:,1), [0 1] ); colormap jet
axis image
yticks( [1  NM/2+1  NM] )
yticklabels( { num2str(NM/2), '0', num2str(NM/2-1) } )
set(gca, 'fontsize', 30, 'fontweight', 'bold' )
ylabel( 'Radius' )
xlabel( 'Projection' )
colorbar

% comparison to Nyquist
DCF_Nyquist = round( (abs((1:NM)-(NM/2+1)).^2)* (4*pi)/2 );
figure(2)
plot( DCF_Nyquist, 'k-' ); 
hold on; plot( sum(gDCF,2), 'r--' ); axis square; grid on
