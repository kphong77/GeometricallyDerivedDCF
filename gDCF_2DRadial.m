
function [DCF] = OptimalDCF_2DRadial( AngleInfo_input, Option ) 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% %
% % Written by KyungPyo Hong
% % Version 1.0:  February, 2019.
% % Version 1.1:  added a weighted-contrast function (KWIC)
% %
% %   - AngleInfo_input:  [0, 2*pi) 
% %       : A x B = projections x repeats (or time frames)
% %   - Option.Nsamples: the number of sampled points in a projection (even number)
% %   - Option.AnglePrecision: round-up for matching up with the precision of AngleInfo_input
% %                            (0 = full precision)
% %   - Option.WeightedContrast: weight-adjust for various image contrasts 
% %                              (Song HK, et al. Magn Reson Med. 2000 Dec;44(6):825-32.) 
% %
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

    % Unit metric
    %  # discrete unit radial/arc distance
    %     - for MRI, dk (Nyquist distance) = 1/FOV = dr (radial) = dL (arc; circumferential) 
    %     - for CT, SPECT, and PET, sensor-to-sensor distance = dr = dL  
    UnitDistance = 1;
    UnitMetric = UnitDistance;
    
    NM = Option.Nsamples;
    [ Nproj  Nrep ] = size( AngleInfo_input );
    
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
    
    if Option.Display == 1
        proj_info = -0.5:1/NM:0.5-1/NM;
        figure(101)
        for ii = 1:Nrep
            for jj = 1:Nproj
                GridInfo = proj_info*exp( 1i*AngleInfo_input(jj,ii) );
                plot( real( GridInfo ), imag( GridInfo ) )
                hold on
            end
            axis square; 
        end
    end

    % prepare angle info
    AngleInfo = mod( [AngleInfo_input;   AngleInfo_input+pi], 2*pi );  
    if CalcPrecision > 1
        AngleInfo = round( AngleInfo * CalcPrecision )/CalcPrecision;
    end
    

    
    DCF = zeros( NM, Nproj, Nrep ); 
    for TimeFrame = 1:Nrep
        fprintf( 'TimeFrame = %3d\n', TimeFrame );
        [tmp_ index_proj ] =  find( WeightedContrast(NM/2+1,:,TimeFrame) > 0 ); 
        NumOfProj = length( index_proj );
        DCF( NM/2+1, index_proj, TimeFrame ) = 1/NumOfProj;   % at radius = 0
        tmp_DCF = zeros( NM/2-1, Nproj*2 );
        parfor radius = 1:NM/2-1
% % %             fprintf( 'radius = %3d\n', radius )
            
            tmp_WC = WeightedContrast( [NM/2+1-radius  NM/2+1+radius], :, TimeFrame )';
            [ index_WC  tmp_ ] = find( tmp_WC(:) > 0 );
            SumOfOverlaps = zeros( 1, Nproj*2 );
            tmp_calc = [];
            for proj_loc = index_WC'
                % % overlap = dL - L
                % % where L = the measured arc distance = radius*|theta2-theta1|
                
                Overlap = [];
                
                ang1 = AngleInfo( proj_loc, TimeFrame );
                ang2 = AngleInfo( index_WC, TimeFrame );
                ang3 = max( ang1, ang2 );
                ang4 = min( ang1, ang2 );
                tmp_calc(:,1) = abs( ang3 - ang4 );
                tmp_calc(:,2) = abs( ang3 - (ang4+2*pi) );
                tmp_sel = sort( tmp_calc, 2, 'ascend' );
                if CalcPrecision > 1
                    tmp_sel = round( tmp_sel*CalcPrecision )/CalcPrecision;
                end
                Overlap = (UnitMetric - radius*tmp_sel(:,1))/UnitMetric;
                [tmp_ind, ~] = find( Overlap <= 0 );
                Overlap( tmp_ind ) = 0;
                SumOfOverlaps(proj_loc) = SumOfOverlaps(proj_loc) + sum(Overlap);
            end
            tmp_DCF( radius, : ) = 1./SumOfOverlaps;
        end
        DCF( NM/2:-1:2, :, TimeFrame ) = tmp_DCF(:,1:Nproj);
        DCF( NM/2+2:end, :, TimeFrame ) = tmp_DCF(:,Nproj+1:end);
    end
    DCF( find( DCF == inf ) ) = 0;
    DCF(1,:,:) = DCF(2,:,:);
    DCF( find( DCF == 0 ) ) = eps;
end
