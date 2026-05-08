#include "mincut_only.h"


int MincutOnly::main() {

    this->WriteToLogFile("Parsing connectedness criterion" , Log::info);
/* F(n) = C log_x(n), where C and x are parameters specified by the user (our default is C=1 and x=10) */
/* G(n) = C n^x, where C and x are parameters specified by the user (here, presumably 0<x<2). Note that x=1 makes it linear. */
        /* static inline bool IsWellConnected(std::string connectedness_criterion, int in_partition_size, int out_partition_size, int edge_cut_size) { */
    double connectedness_criterion_c = this->connectedness_criterion_c;
    double connectedness_criterion_x = this->connectedness_criterion_x;
    double pre_computed_log = this->pre_computed_log;
    ConnectednessCriterion current_connectedness_criterion = this->current_connectedness_criterion;
    std::string connectedness_criterion_custom_string = this->connectedness_criterion_custom_string;
    this->WriteToLogFile("Loading the initial graph" , Log::info);


    std::map<std::string, int> original_to_new_id_map = this->GetOriginalToNewIdMap(this->edgelist);
    std::map<int, std::string> new_to_originial_id_map = this->InvertMap(original_to_new_id_map);
    igraph_t graph;
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);
    this->LoadEdgesFromFile(&graph, this->edgelist, original_to_new_id_map);
    this->WriteToLogFile("Finished loading the initial graph" , Log::info);

    int before_mincut_number_of_clusters = -1;
    int after_mincut_number_of_clusters = -2;
    int iter_count = 0;

    std::map<int, int> new_node_id_to_cluster_id_map = ConstrainedClustering::ReadCommunities(original_to_new_id_map, this->existing_clustering);
    /* std::cerr << "num nodes in mapping: " << original_to_new_id_map.size()  << std::endl; */
    /* std::cerr << "num nodes read from clustering: " << new_node_id_to_cluster_id_map.size()  << std::endl; */
    /* for (auto const& [node_id, cluster_id] : new_node_id_to_cluster_id_map) { */
    /*     std::cerr << std::to_string(node_id) << "," << std::to_string(cluster_id)  << std::endl; */
    /* } */
    ConstrainedClustering::RemoveInterClusterEdges(&graph, new_node_id_to_cluster_id_map);
    /* std::cerr << "num nodes after intercluster removal: " << igraph_vcount(&graph)  << std::endl; */
    /* std::cerr << "num edges after intercluster removal: " << igraph_ecount(&graph)  << std::endl; */

    /** SECTION Get Connected Components START **/
    std::vector<std::vector<int>> connected_components_vector = ConstrainedClustering::GetConnectedComponents(&graph);
    /** SECTION Get Connected Components END **/
    if(current_connectedness_criterion == ConnectednessCriterion::Simple) {
        for(size_t i = 0; i < connected_components_vector.size(); i ++) {
            MincutOnly::done_being_mincut_clusters.push(connected_components_vector[i]);
        }
    } else {
        std::vector<int> top_to_global_map(igraph_vcount(&graph));
        for (int j = 0; j < igraph_vcount(&graph); j++) top_to_global_map[j] = j;
        for (size_t i = 0; i < connected_components_vector.size(); i++) {
            MincutOnly::MinCutWorkerRecursive(
                &graph,
                connected_components_vector[i],
                top_to_global_map,
                current_connectedness_criterion,
                connectedness_criterion_c,
                connectedness_criterion_x,
                pre_computed_log,
                connectedness_criterion_custom_string,
                this->mincut_type);
        }
        this->WriteToLogFile("All top-level CCs processed (recursive driver)", Log::info);
    }


    this->WriteYieldSummary();

    this->WriteToLogFile("Writing output to: " + this->output_file, Log::info);
    this->WriteClusterQueue(MincutOnly::done_being_mincut_clusters, &graph, new_to_originial_id_map);
    MincutOnly::DumpCounters(this);
    igraph_destroy(&graph);
    return 0;
}
