#ifndef CM_H
#define CM_H
#include "constrained.h"
#include <atomic>


class CM : public ConstrainedClustering {
    public:
        CM(std::string edgelist, std::string algorithm, double clustering_parameter, std::string existing_clustering, int num_processors, std::string output_file, std::string log_file, std::string history_file, int log_level, std::string connectedness_criterion, bool prune, std::string mincut_type) : ConstrainedClustering(edgelist, algorithm, clustering_parameter, existing_clustering, num_processors, output_file, log_file, history_file, log_level, connectedness_criterion, mincut_type) {
        };
        int main() override;

        // -------------------- path-trigger counters (per process) --------------------
        // Each path increments its counter from inside the worker(s). Static-inline
        // atomics so all worker threads in this rank share one counter. Dump only at
        // the end of CM::main() via DumpCounters() — matches NJIT-CM's behavior.
        static inline std::atomic<size_t> n_worker_calls{0};       // # times a cluster popped / entered the worker
        static inline std::atomic<size_t> n_mincut_calls{0};       // # invocations of MinCutCustom::ComputeMinCut
        static inline std::atomic<size_t> n_leaf_trims{0};         // # iterations the prune branch peeled a singleton-side
        static inline std::atomic<size_t> n_well_connected{0};     // # clusters that exited as well-connected
        static inline std::atomic<size_t> n_non_trivial_cuts{0};   // # clusters that resulted in a non-trivial mincut
        static inline std::atomic<size_t> n_partitionings{0};      // # leiden re-clustering invocations (RunClusterOnPartition)
        static inline std::atomic<size_t> n_disintegrations{0};    // # clusters that fully disintegrated during peel

        static inline void DumpCounters(ConstrainedClustering* cc) {
            cc->WriteToLogFile("[CM_COUNTERS] worker_calls="     + std::to_string(n_worker_calls.load()),     Log::info);
            cc->WriteToLogFile("[CM_COUNTERS] mincut_calls="     + std::to_string(n_mincut_calls.load()),     Log::info);
            cc->WriteToLogFile("[CM_COUNTERS] leaf_trims="       + std::to_string(n_leaf_trims.load()),       Log::info);
            cc->WriteToLogFile("[CM_COUNTERS] well_connected="   + std::to_string(n_well_connected.load()),   Log::info);
            cc->WriteToLogFile("[CM_COUNTERS] non_trivial_cuts=" + std::to_string(n_non_trivial_cuts.load()), Log::info);
            cc->WriteToLogFile("[CM_COUNTERS] partitionings="    + std::to_string(n_partitionings.load()),    Log::info);
            cc->WriteToLogFile("[CM_COUNTERS] disintegrations="  + std::to_string(n_disintegrations.load()),  Log::info);
        }

        static inline std::vector<std::vector<int>> RunClusterOnPartition(const igraph_t* graph, std::string algorithm, int seed, double clustering_parameter, std::vector<int>& partition) {
            std::vector<std::vector<int>> cluster_vectors;
            std::map<int, std::vector<int>> cluster_map;
            igraph_vector_int_t nodes_to_keep;
            igraph_vector_int_t new_id_to_old_id_map; // here, new_id is the id of the current half, old_id is the id of the induced_subgraph of the cluster (both halves)
            igraph_vector_int_init(&new_id_to_old_id_map, partition.size());
            igraph_vector_int_init(&nodes_to_keep, partition.size());
            for(size_t i = 0; i < partition.size(); i ++) {
                VECTOR(nodes_to_keep)[i] = partition[i];
            }
            igraph_t induced_subgraph;
            igraph_induced_subgraph_map(graph, &induced_subgraph, igraph_vss_vector(&nodes_to_keep), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH, NULL, &new_id_to_old_id_map);
            std::map<int, int> partition_map  = ConstrainedClustering::GetCommunities("", algorithm, seed, clustering_parameter, &induced_subgraph);
            ConstrainedClustering::RemoveInterClusterEdges(&induced_subgraph, partition_map);
            std::vector<std::vector<int>> connected_components_vector = ConstrainedClustering::GetConnectedComponents(&induced_subgraph);
            for(size_t i = 0; i < connected_components_vector.size(); i ++) {
                std::vector<int> translated_cluster_vector;
                for(size_t j = 0; j < connected_components_vector.at(i).size(); j ++) {
                    int new_id = connected_components_vector.at(i).at(j);
                    translated_cluster_vector.push_back(VECTOR(new_id_to_old_id_map)[new_id]);
                }
                cluster_vectors.push_back(translated_cluster_vector);
            }
            igraph_vector_int_destroy(&nodes_to_keep);
            igraph_vector_int_destroy(&new_id_to_old_id_map);
            igraph_destroy(&induced_subgraph);
            return cluster_vectors;
        }

        static inline void MinCutOrClusterWorker(const igraph_t* graph, std::string algorithm, int seed, double clustering_parameter, ConnectednessCriterion current_connectedness_criterion, double connectedness_criterion_c, double connectedness_criterion_x, double pre_computed_log, bool prune, const std::string& connectedness_criterion_custom_string, std::string mincut_type = "cactus") {
            while (true) {
                std::unique_lock<std::mutex> to_be_mincut_lock{to_be_mincut_mutex};
                std::pair<std::vector<int>, int> current_front = CM::to_be_mincut_clusters.front();
                std::vector<int> current_cluster = current_front.first;
                int current_cluster_id = current_front.second;
                std::set<int> current_cluster_set(current_cluster.begin(), current_cluster.end());
                CM::to_be_mincut_clusters.pop();
                to_be_mincut_lock.unlock();
                if(current_cluster.size() == 1 && current_cluster[0] == -1) {
                    // done with work!
                    return;
                }
                CM::n_worker_calls.fetch_add(1, std::memory_order_relaxed);
                igraph_vector_int_t nodes_to_keep;
                igraph_vector_int_t new_id_to_old_id_vector_map;
                igraph_vector_int_init(&nodes_to_keep, current_cluster.size());
                for(size_t i = 0; i < current_cluster.size(); i ++) {
                    VECTOR(nodes_to_keep)[i] = current_cluster[i];
                }
                std::vector<int> in_partition;
                std::vector<int> out_partition;
                bool is_well_connected = false;
                bool is_non_trivial_cut = false;
                int edge_cut_size = -1;
                igraph_t induced_subgraph;
                igraph_vector_int_init(&new_id_to_old_id_vector_map, igraph_vector_int_size(&nodes_to_keep));
                igraph_induced_subgraph_map(graph, &induced_subgraph, igraph_vss_vector(&nodes_to_keep), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH, NULL, &new_id_to_old_id_vector_map);
                std::map<int, int> new_id_to_old_id_map; // new_id = induced subgraph id, old_id = global graph id
                std::map<int, int> old_id_to_new_id_map;
                for(int i = 0; i < igraph_vector_int_size(&new_id_to_old_id_vector_map); i ++) {
                    new_id_to_old_id_map[i] = VECTOR(new_id_to_old_id_vector_map)[i];
                    old_id_to_new_id_map[VECTOR(new_id_to_old_id_vector_map)[i]] = i;
                }
                /* while(!(is_non_trivial_cut || is_well_connected)) { */
                // we do mincuts until it's a non trivial mincut or if the current cluster is already well connected (perhaps after some trivial mincuts)
                while(true) {
                    MinCutCustom mcc(&induced_subgraph, mincut_type);
                    CM::n_mincut_calls.fetch_add(1, std::memory_order_relaxed);
                    edge_cut_size = mcc.ComputeMinCut();
                    std::vector<int> in_partition_candidate = mcc.GetInPartition();
                    std::vector<int> out_partition_candidate = mcc.GetOutPartition();
                    is_well_connected = ConstrainedClustering::IsWellConnected(current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log, in_partition_candidate.size(), out_partition_candidate.size(), edge_cut_size, connectedness_criterion_custom_string);

                    size_t non_trivial_cutoff = 1;
                    is_non_trivial_cut = in_partition_candidate.size() > non_trivial_cutoff && out_partition_candidate.size() > non_trivial_cutoff;
                    if(!prune || is_well_connected || is_non_trivial_cut) {
                        in_partition = in_partition_candidate;
                        out_partition = out_partition_candidate;
                        break;
                    } else {
                        // since we're removing a node, the new_id (induced subgraph id) to old_id mapping needs to change
                        // this is a temp map that will replace the newnew_id_to_old_id_map everytime the induced subgraph is changed
                        std::map<int, int> newnew_id_to_old_id_map;
                        std::vector<int> partitions[] = {in_partition_candidate, out_partition_candidate};
                        int current_cluster_size = in_partition_candidate.size() + out_partition_candidate.size();
                        for(int partition_id = 0; partition_id < 2; partition_id ++) {
                            std::vector<int> current_partition = partitions[partition_id];
                            if(current_partition.size() == 1) {
                                CM::n_leaf_trims.fetch_add(1, std::memory_order_relaxed);
                                int node_to_remove = -1;
                                igraph_vector_int_t newnew_id_to_new_id_map;
                                igraph_vector_int_init(&newnew_id_to_new_id_map, current_cluster_size - 1);
                                node_to_remove = current_partition.at(0);
                                current_cluster_set.erase(new_id_to_old_id_map[node_to_remove]);
                                // MARK: possibly changes behavior?
                                // igraph_delete_vertices_idx(&induced_subgraph, igraph_vss_1(node_to_remove), NULL, &newnew_id_to_new_id_map);
                                igraph_delete_vertices_map(&induced_subgraph, igraph_vss_1(node_to_remove), NULL, &newnew_id_to_new_id_map);
                                for(int i = 0; i < igraph_vector_int_size(&newnew_id_to_new_id_map); i ++) {
                                    int newnew_id = i;
                                    int new_id = VECTOR(newnew_id_to_new_id_map)[newnew_id];
                                    int old_id = new_id_to_old_id_map[new_id];
                                    old_id_to_new_id_map[old_id] = newnew_id;
                                    newnew_id_to_old_id_map[newnew_id] = old_id;
                                }
                                new_id_to_old_id_map = newnew_id_to_old_id_map;
                                igraph_vector_int_destroy(&newnew_id_to_new_id_map);
                                current_cluster_size --;
                            }
                        }
                        if(current_cluster_size == 0) {
                            // this cluster completely disintegrated somehow?
                            CM::n_disintegrations.fetch_add(1, std::memory_order_relaxed);
                            break;
                        }
                    }
                }
                // by this point, the cluster is either well-connected or we actually have a mincut we need to do
                if(is_well_connected) {
                    CM::n_well_connected.fetch_add(1, std::memory_order_relaxed);
                    // just return the cluster without using the mincut
                    current_cluster = std::vector(current_cluster_set.begin(), current_cluster_set.end());
                    {
                        std::lock_guard<std::mutex> done_being_clustered_guard(CM::done_being_clustered_mutex);
                        CM::done_being_clustered_clusters.push({current_cluster, current_cluster_id});
                    }
                } else {
                    CM::n_non_trivial_cuts.fetch_add(1, std::memory_order_relaxed);
                    // do the mincut and actually run a clustering algorithm on both sides
                    std::vector<int> partitions[] = {in_partition, out_partition};
                    for(int partition_id = 0; partition_id < 2; partition_id ++) {
                        std::vector<int> current_partition = partitions[partition_id];
                        // assert(current_partition.size() > 1);
                        if (current_partition.size() > 1) {
                            CM::n_partitionings.fetch_add(1, std::memory_order_relaxed);
                            std::vector<std::vector<int>> current_clusters = CM::RunClusterOnPartition(&induced_subgraph, algorithm, seed, clustering_parameter, current_partition);
                            for(size_t i = 0; i < current_clusters.size(); i ++) {
                                std::vector<int> translated_current_clusters;
                                for(size_t j = 0; j < current_clusters[i].size(); j ++) {
                                    translated_current_clusters.push_back(new_id_to_old_id_map[current_clusters[i][j]]);
                                }
                                {
                                    std::lock_guard<std::mutex> to_be_clustered_guard(CM::to_be_clustered_mutex);
                                    CM::to_be_clustered_clusters.push({translated_current_clusters, current_cluster_id});
                                }
                            }
                        }
                    }
                }
                igraph_vector_int_destroy(&nodes_to_keep);
                igraph_vector_int_destroy(&new_id_to_old_id_vector_map);
                igraph_destroy(&induced_subgraph);
            }
        }

        // ----------------------------------------------------------------------------
        // Recursive variant of MinCutOrClusterWorker.
        //
        // Behaviorally identical to MinCutOrClusterWorker except for two things:
        //   1. Each call's `parent_graph` is the immediate parent's induced subgraph
        //      (NOT the original full graph), so igraph_induced_subgraph_map shrinks
        //      geometrically with depth instead of always rescanning the full input.
        //   2. Sub-cluster recursion happens via direct call (recursion) instead of
        //      through to_be_mincut_clusters / to_be_clustered_clusters queues.
        //
        // Single-threaded only: writes done_being_clustered_clusters and
        // parent_to_child_map without locking. Run with --num-processors 1.
        //
        // Args:
        //   parent_graph              — immediate parent's igraph_t (already restricted)
        //   current_cluster           — node IDs in parent_graph idx space
        //   parent_to_global_map      — parent_graph idx → full-graph idx
        //   current_cluster_id        — id for output / history
        //   next_cluster_id_ref       — counter for assigning sub-cluster ids (mutable)
        //   parent_to_child_map       — written: this cluster's children get registered
        // ----------------------------------------------------------------------------
        static inline void MinCutOrClusterWorkerRecursive(
                const igraph_t* parent_graph,
                const std::vector<int>& current_cluster,
                const std::vector<int>& parent_to_global_map,
                int current_cluster_id,
                int& next_cluster_id_ref,
                std::map<int, std::vector<int>>* parent_to_child_map,
                std::string algorithm,
                int seed,
                double clustering_parameter,
                ConnectednessCriterion current_connectedness_criterion,
                double connectedness_criterion_c,
                double connectedness_criterion_x,
                double pre_computed_log,
                bool prune,
                const std::string& connectedness_criterion_custom_string,
                std::string mincut_type = "cactus") {
            CM::n_worker_calls.fetch_add(1, std::memory_order_relaxed);

            if (current_cluster.size() < 2) {
                std::vector<int> done_cluster;
                done_cluster.reserve(current_cluster.size());
                for (int idx : current_cluster) done_cluster.push_back(parent_to_global_map[idx]);
                CM::done_being_clustered_clusters.push({done_cluster, current_cluster_id});
                CM::n_well_connected.fetch_add(1, std::memory_order_relaxed);
                return;
            }

            // Build induced subgraph from parent_graph using current_cluster.
            igraph_vector_int_t nodes_to_keep;
            igraph_vector_int_init(&nodes_to_keep, current_cluster.size());
            for (size_t i = 0; i < current_cluster.size(); i++) {
                VECTOR(nodes_to_keep)[i] = current_cluster[i];
            }
            igraph_vector_int_t new_id_to_parent_id_vec;
            igraph_vector_int_init(&new_id_to_parent_id_vec, current_cluster.size());
            igraph_t induced_subgraph;
            igraph_induced_subgraph_map(parent_graph, &induced_subgraph,
                igraph_vss_vector(&nodes_to_keep),
                IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH,
                NULL, &new_id_to_parent_id_vec);

            // live_to_global: induced_subgraph idx → full-graph idx. Updated as nodes
            // are peeled by singleton-side trims.
            std::vector<int> live_to_global(igraph_vcount(&induced_subgraph));
            for (int i = 0; i < igraph_vcount(&induced_subgraph); i++) {
                int parent_idx = VECTOR(new_id_to_parent_id_vec)[i];
                live_to_global[i] = parent_to_global_map[parent_idx];
            }

            int current_cluster_size = igraph_vcount(&induced_subgraph);
            std::vector<int> in_partition, out_partition;
            bool is_well_connected = false;
            bool is_non_trivial_cut = false;
            int edge_cut_size = -1;

            // Mincut + (prune=true only) singleton-side trim loop.
            // Break condition mirrors the original worker: prune=false → take the
            // first cut as-is. prune=true → keep peeling until well-connected or
            // non-trivial.
            while (true) {
                MinCutCustom mcc(&induced_subgraph, mincut_type);
                CM::n_mincut_calls.fetch_add(1, std::memory_order_relaxed);
                edge_cut_size = mcc.ComputeMinCut();
                std::vector<int> in_partition_candidate  = mcc.GetInPartition();
                std::vector<int> out_partition_candidate = mcc.GetOutPartition();
                is_well_connected = ConstrainedClustering::IsWellConnected(
                    current_connectedness_criterion,
                    connectedness_criterion_c, connectedness_criterion_x,
                    pre_computed_log,
                    in_partition_candidate.size(), out_partition_candidate.size(),
                    edge_cut_size,
                    connectedness_criterion_custom_string);
                is_non_trivial_cut = in_partition_candidate.size() > 1
                                    && out_partition_candidate.size() > 1;

                if (!prune || is_well_connected || is_non_trivial_cut) {
                    in_partition  = in_partition_candidate;
                    out_partition = out_partition_candidate;
                    break;
                }

                // Singleton-side trim: peel the leaf node, remap live_to_global.
                std::vector<int> partitions[] = {in_partition_candidate, out_partition_candidate};
                for (int pid = 0; pid < 2; pid++) {
                    if (partitions[pid].size() == 1) {
                        CM::n_leaf_trims.fetch_add(1, std::memory_order_relaxed);
                        int node_to_remove = partitions[pid].at(0);
                        igraph_vector_int_t newnew_to_new;
                        igraph_vector_int_init(&newnew_to_new, current_cluster_size - 1);
                        igraph_delete_vertices_map(&induced_subgraph, igraph_vss_1(node_to_remove),
                            NULL, &newnew_to_new);
                        std::vector<int> new_live_to_global(current_cluster_size - 1);
                        for (int i = 0; i < current_cluster_size - 1; i++) {
                            int old_idx = VECTOR(newnew_to_new)[i];
                            new_live_to_global[i] = live_to_global[old_idx];
                        }
                        live_to_global = std::move(new_live_to_global);
                        igraph_vector_int_destroy(&newnew_to_new);
                        current_cluster_size--;
                    }
                }
                if (current_cluster_size == 0) {
                    CM::n_disintegrations.fetch_add(1, std::memory_order_relaxed);
                    break;
                }
            }

            if (is_well_connected) {
                CM::n_well_connected.fetch_add(1, std::memory_order_relaxed);
                std::vector<int> done_cluster(live_to_global.begin(), live_to_global.end());
                CM::done_being_clustered_clusters.push({done_cluster, current_cluster_id});
            } else {
                CM::n_non_trivial_cuts.fetch_add(1, std::memory_order_relaxed);
                std::vector<int> partitions[] = {in_partition, out_partition};
                for (int pid = 0; pid < 2; pid++) {
                    const std::vector<int>& cur_part = partitions[pid];
                    if (cur_part.size() <= 1) continue;
                    CM::n_partitionings.fetch_add(1, std::memory_order_relaxed);
                    std::vector<std::vector<int>> ccs = CM::RunClusterOnPartition(
                        &induced_subgraph, algorithm, seed, clustering_parameter,
                        const_cast<std::vector<int>&>(cur_part));
                    for (auto& cc : ccs) {
                        int new_id = next_cluster_id_ref++;
                        (*parent_to_child_map)[current_cluster_id].push_back(new_id);
                        // Recurse: the new parent is THIS induced_subgraph; live_to_global
                        // composes the chain from local idx back to full-graph idx.
                        CM::MinCutOrClusterWorkerRecursive(
                            &induced_subgraph,
                            cc,
                            live_to_global,
                            new_id,
                            next_cluster_id_ref,
                            parent_to_child_map,
                            algorithm, seed, clustering_parameter,
                            current_connectedness_criterion,
                            connectedness_criterion_c, connectedness_criterion_x,
                            pre_computed_log,
                            prune,
                            connectedness_criterion_custom_string,
                            mincut_type);
                    }
                }
            }

            igraph_vector_int_destroy(&nodes_to_keep);
            igraph_vector_int_destroy(&new_id_to_parent_id_vec);
            igraph_destroy(&induced_subgraph);
        }

    private:
        bool prune;
        static inline std::mutex to_be_mincut_mutex;
        static inline std::condition_variable to_be_mincut_condition_variable;
        static inline std::queue<std::pair<std::vector<int>, int>> to_be_mincut_clusters;
        static inline std::mutex to_be_clustered_mutex;
        static inline std::queue<std::pair<std::vector<int>, int>> to_be_clustered_clusters;
        static inline std::mutex done_being_clustered_mutex;
        static inline std::queue<std::pair<std::vector<int>, int>> done_being_clustered_clusters;
};

#endif
