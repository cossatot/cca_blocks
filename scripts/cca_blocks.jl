using Revise

using Oiler

using CSV
using JSON
using DataFrames:eachrow
using DataFrames, DataFramesMeta
using Setfield

using PyPlot

# options
geol_slip_rate_weight = 2.
save_results = true

# load data

# cca
cca_block_file = "../block_data/cca_blocks.geojson"
cca_fault_file = "../block_data/cca_faults.geojson"
cca_slip_rates_file = "../block_data/cca_geol_slip_rates.geojson"

ant_tris_file = "../../subduction/sub_tri_meshes/ant_slab2.geojson"
cam_tris_file = "../../subduction/sub_tri_meshes/cam_slab2_fine.geojson"

cca_bounds_file = "../block_data/cca_block_bounds.geojson"
cca_nsam_bounds_file = "../block_data/cca_nsam_block_bounds.geojson"

# nsam
sam_block_file ="../../sam_blocks/block_data/sam_blocks.geojson"
sam_fault_file ="../../sam_blocks/block_data/sam_faults.geojson"
nsam_tris_file = "../../sam_blocks/block_data/nsam_tris.geojson"
mora_vels_file = "../../sam_blocks/block_data/mora_vels.geojson"

# glo
glo_block_file = "../../global_scale_plates/global_scale_plates.geojson"
glo_fault_file = "../../global_scale_plates/global_scale_faults.geojson"
glo_slip_rates_file = "../../global_scale_plates/global_scale_slip_rates.geojson"
gsrm_vels_file = "../../../gsrm/gps/gps_na.geojson"
midas_vels_file = "../../cascadia_blocks/data/midas_vels.geojson"
garnier_vels_file = "../geod_data/garnier_et_al_2022_vels_igs08.geojson"


@info "joining blocks"
cca_blocks = Oiler.IO.gis_vec_file_to_df(cca_block_file; fid_drop=[])
glo_blocks = Oiler.IO.gis_vec_file_to_df(glo_block_file; fid_drop=["ant"])
sam_blocks = Oiler.IO.gis_vec_file_to_df(sam_block_file)

block_df = vcat(cca_blocks,
                sam_blocks,
                glo_blocks;
                cols=:union)

#@info "removing antarctica for now"
#ant_df = filter(row -> (row.fid == "ant"), block_df)
#block_df = filter(row -> !(row.fid == "ant"), block_df)

println("n blocks: ", size(block_df, 1))

@info "culling blocks"
println("n blocks before ", size(block_df, 1))
bound_df = Oiler.IO.gis_vec_file_to_df(cca_nsam_bounds_file)
block_df = Oiler.IO.get_blocks_in_bounds!(block_df, bound_df; epsg=3995)
println("n blocks after ", size(block_df, 1))

@info "doing faults"
fault_df, faults, fault_vels = Oiler.IO.process_faults_from_gis_files(
                                            cca_fault_file,
                                            sam_fault_file,
                                            glo_fault_file;
                                            block_df=block_df,
                                            subset_in_bounds=true,
                                            #usd_default=1.,
                                            #lsd_default=4.,
                                            e_default=5.,
                                            #fid_drop="ccaf002",
                                            )

println("n faults: ", length(faults))
println("n fault vels: ", length(fault_vels))


@info "doing geologic slip rates"
glo_slip_rate_df = Oiler.IO.gis_vec_file_to_df(glo_slip_rates_file)
cca_slip_rate_df = Oiler.IO.gis_vec_file_to_df(cca_slip_rates_file)

geol_slip_rate_df = vcat(glo_slip_rate_df,
                         cca_slip_rate_df,
                         )

geol_slip_rate_df, geol_slip_rate_vels = Oiler.IO.make_geol_slip_rate_vels!(
                                            geol_slip_rate_df,
                                            fault_df;
                                            weight=geol_slip_rate_weight
                                            )

println("n geol slip rates: ", length(geol_slip_rate_vels))

@info "doing GNSS"
gsrm_vel_df = Oiler.IO.gis_vec_file_to_df(gsrm_vels_file)
midas_vel_df = Oiler.IO.gis_vec_file_to_df(midas_vels_file)
garn_vel_df = Oiler.IO.gis_vec_file_to_df(garnier_vels_file)
mora_vel_df = Oiler.IO.gis_vec_file_to_df(mora_vels_file)

gnss_vel_df = vcat(gsrm_vel_df, midas_vel_df, cols=:union)

@time gsmd_vels = Oiler.IO.make_vels_from_gnss_and_blocks(gnss_vel_df, block_df;
    fix="1111", epsg=102016,
    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err, name=:station
)

@time garn_vels = Oiler.IO.make_vels_from_gnss_and_blocks(garn_vel_df, block_df;
    fix="igs08", epsg=102016,
    ve=:e_vel, vn=:n_vel, ee=:e_rr, en=:n_err, name=:site
)

@time mora_vels = Oiler.IO.make_vels_from_gnss_and_blocks(mora_vel_df, block_df;
    fix="itrf14", epsg=102016,
    ve=:e_vel, vn=:n_vel, ee=:e_err, en=:n_err, name=:station
)
gnss_vels = vcat(gsmd_vels, 
                 garn_vels, 
                 #mora_vels,
                 )

println("n gnss vels: ", length(gnss_vels))

@info "doing tris"
ant_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(ant_tris_file))
cam_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(cam_tris_file))
nsam_tris = Oiler.IO.tris_from_geojson(JSON.parsefile(nsam_tris_file))

tris = vcat(cam_tris,
            #ant_tris,
            nsam_tris
            )

#tris = cam_tris
#tris = []

println("n tris: ", length(tris))

vels = vcat(fault_vels,
            gnss_vels,
            geol_slip_rate_vels
            )

println("n total vels: ", length(vels))
vel_groups = Oiler.group_vels_by_fix_mov(vels)

tri_distance_weight = 5.


@info "Solving"
@time results = Oiler.solve_block_invs_from_vel_groups(vel_groups,
            tris=tris,
            faults=faults,
            elastic_floor=1e-5,
            tri_distance_weight=tri_distance_weight,
            regularize_tris=true,
            tri_priors=false,
            predict_vels=true,
            pred_se=true,
            check_closures=false,
            check_nans=true,
            sparse_lhs=true,
            constraint_method="kkt_sym",
            factorization="lu")


Oiler.WebViewer.write_web_viewer(results=results, block_df=block_df,
                                 directory="../web_viewer", ref_pole="1111")

Oiler.ResultsAnalysis.compare_data_results(results=results,
                                           vel_groups=vel_groups,
                                           geol_slip_rate_df=geol_slip_rate_df,
                                           geol_slip_rate_vels=geol_slip_rate_vels,
                                           fault_df=fault_df)

map_fig = Oiler.Plots.plot_results_map(results, vel_groups, faults, tris)
rates_fig = Oiler.Plots.plot_slip_rate_fig(geol_slip_rate_df, 
                                           geol_slip_rate_vels, 
                                           fault_df, results)

show()

if save_results == true
    Oiler.IO.write_tri_results_to_gj(tris, results,
                                     "../results/cca_tris.geojson",
                                     name="cca tri results")
    Oiler.IO.write_fault_results_to_gj(results,
                                       "../results/cca_nsam_faults.geojson",
                                       name="cca_nsam_faults")
    Oiler.IO.write_gnss_vel_results_to_csv(results, vel_groups;
                                       name="../results/cca_nsam_gnss_results.csv")
end

println("done!")
