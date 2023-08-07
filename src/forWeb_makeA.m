function [A, E_geo, k_mat] = forWeb_makeA(k_factors, sim_type, Lriver_FHgP, IHgD_pristine, IHgP_pristine)
    %==========================================================================
    % OBJECTIVE
    %   Assemble matrix A, such that dM/dt = A*M + E
    %
    % REVISION HISTORY
    %   18 Jan 2012 - hma - version 1.0 written
    %   19 Jul 2012 - hma - remove deep mineral pool and treat geogenic
    %                       emissions as an external forcing
    %   20 Mar 2014 - hma - update parameterization for rivers and add burial
    %                       in coastal benthic sediment
    %   08 Sep 2014 - HMA - clean up code and comments for public release
    %   08 Sep 2014 - HMA - clean up code and comments for public release
    %
    % Helen M. Amos, hamos@hsph.harvard.edu
    %==========================================================================

    % clear variable, for safety's sake
    clear A;

    % First-order rate coefficients, k (1/yr)
    k_mat = forWeb_rate_coeffs(k_factors, Lriver_FHgP, IHgD_pristine, IHgP_pristine);
    
    % translate k values to correct coefficients
    k_A_oHgII = k_mat(1);
    k_A_tHgII = k_mat(2);
    k_A_tHg0 = k_mat(3);
    k_A_oHg0 = k_mat(4);
    fdep_tf = k_mat(5);
    fdep_ts = k_mat(6);
    fdep_ta = k_mat(7);
    k_Oc_ev = k_mat(8);
    k_Oc_sp1 = k_mat(9);
    k_Oc_vsi = k_mat(10);
    k_Oc_sp2 = k_mat(11);
    k_Oc_vis = k_mat(12);
    k_Oc_vid = k_mat(13);
    k_Oc_sp3 = k_mat(14);
    k_Oc_vdi = k_mat(15);
    k_Te_rf = k_mat(16);
    k_Te_p = k_mat(17);
    k_T_exfs = k_mat(18);
    k_T_exfa = k_mat(19);
    k_Te_BBf_1 = k_mat(20);
    k_Te_BBf_2 = k_mat(21);
    k_Te_rs = k_mat(22);
    k_T_exsf = k_mat(23);
    k_T_exsa = k_mat(24);
    k_Te_BBs_1 = k_mat(25);
    k_Te_BBs_2 = k_mat(26);
    k_Te_ra = k_mat(27);
    k_T_exaf = k_mat(28);
    k_T_exam = k_mat(29);
    k_Te_BBa_1 = k_mat(30);
    k_Te_BBa_2 = k_mat(31);
    k_T_riv_f = k_mat(32);
    k_T_riv_s = k_mat(33);
    k_T_riv_a = k_mat(34);
    k_O_riv_f = k_mat(35);
    k_O_riv_s = k_mat(36);
    k_O_riv_a = k_mat(37);
    E_geo = k_mat(38);
    f_HgPexport = k_mat(39);

    % The rate coefficient for biomass burning differs between the
    % anthropogenic and pre-anthropogenic simulations. It's the only rate
    % coefficient with a time-dependence. 
    if sim_type == 1;     % <-- pre-anthropogenic era
        k_Te_BBf = k_Te_BBf_1; % fast soil
        k_Te_BBs = k_Te_BBs_1; % slow soil
        k_Te_BBa = k_Te_BBa_1; % armored soil

    elseif sim_type == 2; % <-- anthropogenic era
        k_Te_BBf = k_Te_BBf_2; % fast soil
        k_Te_BBs = k_Te_BBs_2; % slow soil
        k_Te_BBa = k_Te_BBa_2; % armored soil    
    else
        message('Invalid simulation type! Must be 1 or 2.')
    end

    %--------------------------------------------------------------------------
    % Assemble matrix A, such that dM/dt = A*M + E
    %--------------------------------------------------------------------------

    %-- atmosphere
    Aatm = [-(k_A_oHgII + k_A_oHg0 + k_A_tHgII + k_A_tHg0)   % Matm term
             (k_Te_rf + k_Te_p + k_Te_BBf)                   % Mtf term
             (k_Te_rs + k_Te_BBs)                            % Mts
             (k_Te_ra + k_Te_BBa)                            % Mta
              k_Oc_ev                                        % Mocs
              0                                              % Moci
              0 ];                                           % Mocd

    %-- fast soil pool
    Atf = [ (k_A_tHgII * fdep_tf + k_A_tHg0)
           -(k_T_riv_f + k_Te_rf + k_Te_p + k_T_exfs + k_T_exfa + k_Te_BBf);
             k_T_exsf
             k_T_exaf
             0
             0
             0 ];

    %-- slow soil pool
    Ats = [ k_A_tHgII*fdep_ts
            k_T_exfs 
          -(k_Te_rs + k_T_exsf + k_T_exsa + k_T_riv_s + k_Te_BBs)
            0
            0
            0
            0 ];

    %-- amored soil pool
    Ata = [ k_A_tHgII * fdep_ta
            k_T_exfa
            k_T_exsa
          -(k_Te_ra + k_T_exaf + k_T_riv_a + k_Te_BBa)
            0
            0
            0 ];

    %-- surface ocean
    Aocs = [ k_A_oHgII + k_A_oHg0
             k_O_riv_f
             k_O_riv_s
             k_O_riv_a
           -(k_Oc_sp1 + k_Oc_ev + k_Oc_vsi)
             k_Oc_vis
             0];

    %-- subsurface ocean
    Aoci = [ 0
             0
             0
             0
             k_Oc_sp1 + k_Oc_vsi
           -(k_Oc_vis+ k_Oc_vid + k_Oc_sp2)
             k_Oc_vdi ];

    %-- deep ocean
    Aocd = [ 0
             0
             0
             0
             0
             k_Oc_vid + k_Oc_sp2
           -(k_Oc_vdi + k_Oc_sp3)];


    %-- matrix A     
    A = [Aatm.' ; Atf.' ; Ats.' ; Ata.' ; Aocs.' ; Aoci.' ; Aocd.'];
end