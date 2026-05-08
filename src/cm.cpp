#include "cm.h"

int CM::main() {
    this->WriteToLogFile("Loading the initial graph" , Log::info);
    std::map<std::string, int> original_to_new_id_map = this->GetOriginalToNewIdMap(this->edgelist);
    std::map<int, std::string> new_to_originial_id_map = this->InvertMap(original_to_new_id_map);
    igraph_t graph;
    /* igraph_set_attribute_table(&igraph_cattribute_table); */
    igraph_empty(&graph, 0, IGRAPH_UNDIRECTED);

    this->LoadEdgesFromFile(&graph, this->edgelist, original_to_new_id_map);
    this->WriteToLogFile("Finished loading the initial graph" , Log::info);

    int after_mincut_number_of_clusters = -2;
    int iter_count = 0;
    
    std::map<int, int> node_id_to_cluster_id_map;
    if(this->existing_clustering == "") {
        int seed = 0;
        node_id_to_cluster_id_map = ConstrainedClustering::GetCommunities("", this->algorithm, seed, this->clustering_parameter, &graph);
        ConstrainedClustering::RemoveInterClusterEdges(&graph, node_id_to_cluster_id_map);
    } else if(this->existing_clustering != "") {
        node_id_to_cluster_id_map = ConstrainedClustering::ReadCommunities(original_to_new_id_map, this->existing_clustering);
        ConstrainedClustering::RemoveInterClusterEdges(&graph, node_id_to_cluster_id_map);
    }
    std::map<int, std::vector<int>> cluster_id_to_node_id_map = ConstrainedClustering::ReverseMap(node_id_to_cluster_id_map); // graph node id
    int next_cluster_id = ConstrainedClustering::FindMaxClusterId(cluster_id_to_node_id_map) + 1;

    std::map<int, std::vector<int>> parent_to_child_map;
    std::vector<std::vector<int>> connected_components_vector = ConstrainedClustering::GetConnectedComponents(&graph);
    for(size_t i = 0; i < connected_components_vector.size(); i ++) {
        int parent_cluster_id = -1; // -1 indicates root
        int current_cluster_id = -1;

        int first_node_in_subcluster = connected_components_vector[i][0];
        int original_cluster_id_of_first_node_in_subcluster = node_id_to_cluster_id_map.at(first_node_in_subcluster);
        int size_of_original_cluster = cluster_id_to_node_id_map[original_cluster_id_of_first_node_in_subcluster].size();
        int size_of_sub_cluster = connected_components_vector[i].size();
        if (size_of_original_cluster == size_of_sub_cluster) {
            current_cluster_id = original_cluster_id_of_first_node_in_subcluster;
        } else {
            bool original_cluster_id_of_first_node_in_subcluster_found = false;
            for (size_t j = 0; j < parent_to_child_map[-1].size(); j ++) {
                if (parent_to_child_map[-1][j] == original_cluster_id_of_first_node_in_subcluster) {
                    original_cluster_id_of_first_node_in_subcluster_found = true;
                    break;
                }
            }
            if (!original_cluster_id_of_first_node_in_subcluster_found) {
                parent_to_child_map[-1].push_back(original_cluster_id_of_first_node_in_subcluster);
            }
            parent_cluster_id = original_cluster_id_of_first_node_in_subcluster;
            current_cluster_id = next_cluster_id;
            next_cluster_id ++;
        }
        parent_to_child_map[parent_cluster_id].push_back(current_cluster_id);
        // Recursive driver instead of the queue-based outer loop. Each top-level CC
        // is processed synchronously; each recursive call's parent_graph is the
        // immediate parent's induced subgraph (NOT the original full graph), so
        // igraph_induced_subgraph_map shrinks geometrically with depth.
        // Single-threaded: must run with --num-processors 1.
        // Top-level parent_to_global_map is the identity over the full graph.
        std::vector<int> top_to_global_map(igraph_vcount(&graph));
        for (int j = 0; j < igraph_vcount(&graph); j++) top_to_global_map[j] = j;
        int seed = 0;
        CM::MinCutOrClusterWorkerRecursive(
            &graph,
            connected_components_vector[i],
            top_to_global_map,
            current_cluster_id,
            next_cluster_id,
            &parent_to_child_map,
            this->algorithm,
            seed,
            this->clustering_parameter,
            this->current_connectedness_criterion,
            this->connectedness_criterion_c,
            this->connectedness_criterion_x,
            this->pre_computed_log,
            this->prune,
            this->connectedness_criterion_custom_string,
            this->mincut_type);
    }
    this->WriteToLogFile("All top-level CCs processed (recursive driver)", Log::info);
    iter_count = 1;
    after_mincut_number_of_clusters = 0;


    this->WriteYieldSummary();

    this->WriteToLogFile("Writing output to: " + this->output_file, Log::info);
    this->WriteClusterQueue(CM::done_being_clustered_clusters, &graph, new_to_originial_id_map);
    this->WriteClusterHistory(parent_to_child_map);
    CM::DumpCounters(this);
    igraph_destroy(&graph);
    return 0;
}
