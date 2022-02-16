
function [DCF_3D] = gDCF_3DRadial( AngleInfo_input, Option ) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Written by KyungPyo Hong 
% % Version 1.0:  March, 2019.  based on surface area
% % Version 1.1:  added a weighted-contrast function (KWIC)
% %
% %  - AngleInfo_input (radian), 3-dimension = A x B x C
% %      A: projections
% %      B: 1=azimuthal angle [0,2pi); 2=polar angle [0,pi)
% %      C: repeats (or time frames) 
% %  - Option.Nsamples: the number of sampled points in a projection (even number)
% %  - Option.AnglePrecision: round-up for matching up with the precision of AngleInfo_input
% %                           (0 = full precision)
% %  - Option.WeightedContrast: weight-adjust for various image contrasts 
% %                             (Song HK, et al. Magn Reson Med. 2000 Dec;44(6):825-32.) 
% %
% % point R = (radius r, azimuthal theta, polar phi)
% %  --> P(x,y,z):  x = r*sin(phi)*cos(theta)
% %                 y = r*sin(phi)*sin(theta)
% %                 z = r*cos(phi)
% % *calculation of an angular distance between two points on the sphere
% %    - Pa = (x1,y1,z1)
% %    - Pb = (x2,y2,z2)
% %    - angular distance = atan2( norm(cross(Pa,Pb)), dot(Pa,Pb) )
% %
% % Publication: "Accelerating compressed sensing reconstruction of subsampled radial k-space data using geometrically-derived density compensation"
% % Hong K, et al. Physics in Medicine & Biology, Volume 66, Number 21.
% %  - Citation: KyungPyo Hong et al 2021 Phys. Med. Biol. 66 21NT01
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Unit metric
    %   # for 3D radial (surface area)
    %       - Unit metric = unit area = (unit distance)^2
    UnitDistance = 1;
    UnitMetric = UnitDistance^2;
    
    NM = Option.Nsamples;
    [ Nproj  Nangles  Nrep ] = size( AngleInfo_input );
    
    if ~isfield( Option, 'AnglePrecision' ) 
        Option.AnglePrecision = 0;
    end
    CalcPrecision = 10^(Option.AnglePrecision);
    
    if ~isfield( Option, 'WeightedContrast' ) || length( Option.WeightedContrast(:) ) == 1
        WeightedContrast = ones( NM, Nproj, Nrep );
    else
        WeightedContrast = Option.WeightedContrast;
    end
    
    if ~isfield( Option, 'Display' )
        Option.Display = 0;
    end
    
    proj_info = -0.5:1/NM:0.5-1/NM;
    proj_grid = zeros( NM, Nproj, Nrep, 3 );
    for ii = 1:Nrep
        proj_grid(:,:,ii,1) = (proj_info.*cos(AngleInfo_input(:,1,ii)).*sin(AngleInfo_input(:,2,ii)))'; 
        proj_grid(:,:,ii,2) = (proj_info.*sin(AngleInfo_input(:,1,ii)).*sin(AngleInfo_input(:,2,ii)))'; 
        proj_grid(:,:,ii,3) = (proj_info.*cos(AngleInfo_input(:,2,ii)))'; 
    end
    proj_grid = proj_grid*(-1);
    if Option.Display == 1
        figure(101)
        for ii = 1:Nrep
            for jj = 1:Nproj
                plot3( proj_grid(:,jj,ii,1), proj_grid(:,jj,ii,2), proj_grid(:,jj,ii,3) ) 
                hold on
            end
        end
        grid on
    end
    
   
    DCF_3D = zeros( NM, Nproj, Nrep );
    for TimeFrame = 1:Nrep
        fprintf( 'TimeFrame = %3d\n', TimeFrame );
        [tmp_ index_proj ] =  find( WeightedContrast(NM/2+1,:,TimeFrame) > 0 ); 
        NumOfProj = length( index_proj );
        DCF_3D( NM/2+1, index_proj, TimeFrame ) = 1/NumOfProj;   % at radius = 0

        tmp_DCF_3D = zeros( NM/2-1, Nproj*2 );
        parfor radius = 1:NM/2-1
            fprintf( 'radius = %3d\n', radius )
            VectorPos = squeeze(proj_grid(NM/2+1+radius,:,TimeFrame,:));
            VectorNeg = squeeze(proj_grid(NM/2+1-radius,:,TimeFrame,:));
            AllVectors = [VectorPos; VectorNeg];
            
            tmp_WC = WeightedContrast( [NM/2+1-radius  NM/2+1+radius], :, TimeFrame )';
            [ index_WC  tmp_ ] = find( tmp_WC(:) > 0 );

            SumOfOverlaps = zeros( 1, Nproj*2 );
            for proj_loc = index_WC'
                % % overlap for 3D radial: based on the surface area of sphere 
                % %  - dA = Unit Metric = unit area = (unit distance)^2
                % %  -  A = radius^2 * (theta2-theta1) * |cos(phi1)-cos(phi2)|
                % %       = 2*pi*(radius^2)*(1-cos(AngularDistanceBtwTwoPoints)) at the pole (r,0,0)
                
                % % calculate an angular distance between two points on the surface at a given radius
                dotAB = sum( AllVectors(proj_loc,:).*AllVectors(index_WC,:), 2 );
                tmp_AngDist = atan2( sqrt( sum( cross( ones(length(index_WC),3).*AllVectors(proj_loc,:), AllVectors(index_WC,:) ).^2,2)), dotAB );
                AngularDistance = abs( tmp_AngDist );
                if CalcPrecision > 1
                    AngularDistance = round( AngularDistance * CalcPrecision )/CalcPrecision;
                end

                Overlap = (UnitMetric - 2*pi*(radius^2)*(1-cos(AngularDistance)))/UnitMetric;
                [tmp_ind, ~] = find( Overlap <= 0 );
                Overlap( tmp_ind ) = 0;
                SumOfOverlaps(proj_loc) = SumOfOverlaps(proj_loc) + sum(Overlap);
            end
            tmp_DCF_3D( radius, : ) = 1./SumOfOverlaps;
        end
        DCF_3D(NM/2:-1:2,:,TimeFrame) = tmp_DCF_3D(:,1:Nproj);
        DCF_3D(NM/2+2:end,:,TimeFrame) = tmp_DCF_3D(:,Nproj+1:end);
    end
    DCF_3D( find( DCF_3D == inf ) ) = 0;
    DCF_3D(1,:,:) = DCF_3D(2,:,:);
    DCF_3D( find( DCF_3D == 0 ) ) = eps;
end


