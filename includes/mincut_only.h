#ifndef MINCUT_ONLY_H
#define MINCUT_ONLY_H
#include "constrained.h"
#include <atomic>


class MincutOnly : public ConstrainedClustering {
    public:
        MincutOnly(std::string edgelist, std::string existing_clustering, int num_processors, std::string output_file, std::string log_file, int log_level, std::string connectedness_criterion, std::string mincut_type) : ConstrainedClustering(edgelist, "", -1, existing_clustering, num_processors, output_file, log_file, "", log_level, connectedness_criterion, mincut_type) {
        };
        int main() override;

        static inline std::atomic<size_t> n_worker_calls{0};
        static inline std::atomic<size_t> n_mincut_calls{0};
        static inline std::atomic<size_t> n_well_connected{0};
        static inline std::atomic<size_t> n_non_trivial_cuts{0};
        static inline std::atomic<size_t> n_partitionings{0};

        static inline void DumpCounters(ConstrainedClustering* cc) {
            cc->WriteToLogFile("[MO_COUNTERS] worker_calls="     + std::to_string(n_worker_calls.load()),     Log::info);
            cc->WriteToLogFile("[MO_COUNTERS] mincut_calls="     + std::to_string(n_mincut_calls.load()),     Log::info);
            cc->WriteToLogFile("[MO_COUNTERS] well_connected="   + std::to_string(n_well_connected.load()),   Log::info);
            cc->WriteToLogFile("[MO_COUNTERS] non_trivial_cuts=" + std::to_string(n_non_trivial_cuts.load()), Log::info);
            cc->WriteToLogFile("[MO_COUNTERS] partitionings="    + std::to_string(n_partitionings.load()),    Log::info);
        }

        static inline std::vector<std::vector<int>> GetConnectedComponentsOnPartition(const igraph_t* graph, std::vector<int>& partition) {
            std::vector<std::vector<int>> cluster_vectors;
            igraph_vector_int_t nodes_to_keep;
            igraph_vector_int_t new_id_to_old_id_map;
            igraph_vector_int_init(&new_id_to_old_id_map, partition.size());
            igraph_vector_int_init(&nodes_to_keep, partition.size());
            for(size_t i = 0; i < partition.size(); i ++) {
                VECTOR(nodes_to_keep)[i] = partition[i];
            }
            igraph_t induced_subgraph;
            igraph_induced_subgraph_map(graph, &induced_subgraph, igraph_vss_vector(&nodes_to_keep), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH, NULL, &new_id_to_old_id_map);
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

        static inline void MinCutWorkerRecursive(
                const igraph_t* parent_graph,
                const std::vector<int>& current_cluster,
                const std::vector<int>& parent_to_global_map,
                ConnectednessCriterion current_connectedness_criterion,
                double connectedness_criterion_c,
                double connectedness_criterion_x,
                double pre_computed_log,
                const std::string& connectedness_criterion_custom_string,
                std::string mincut_type = "cactus") {
            MincutOnly::n_worker_calls.fetch_add(1, std::memory_order_relaxed);

            if (current_cluster.size() < 2) {
                std::vector<int> done_cluster;
                done_cluster.reserve(current_cluster.size());
                for (int idx : current_cluster) done_cluster.push_back(parent_to_global_map[idx]);
                MincutOnly::done_being_mincut_clusters.push(done_cluster);
                MincutOnly::n_well_connected.fetch_add(1, std::memory_order_relaxed);
                return;
            }

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

            std::vector<int> child_to_global(igraph_vcount(&induced_subgraph));
            for (int i = 0; i < igraph_vcount(&induced_subgraph); i++) {
                int parent_idx = VECTOR(new_id_to_parent_id_vec)[i];
                child_to_global[i] = parent_to_global_map[parent_idx];
            }

            MinCutCustom mcc(&induced_subgraph, mincut_type);
            MincutOnly::n_mincut_calls.fetch_add(1, std::memory_order_relaxed);
            int edge_cut_size = mcc.ComputeMinCut();
            std::vector<int> in_partition  = mcc.GetInPartition();
            std::vector<int> out_partition = mcc.GetOutPartition();
            bool current_criterion = ConstrainedClustering::IsWellConnected(
                current_connectedness_criterion,
                connectedness_criterion_c, connectedness_criterion_x,
                pre_computed_log,
                in_partition.size(), out_partition.size(), edge_cut_size,
                connectedness_criterion_custom_string);

            if (!current_criterion) {
                MincutOnly::n_non_trivial_cuts.fetch_add(1, std::memory_order_relaxed);
                std::vector<int> partitions[] = {in_partition, out_partition};
                for (int pid = 0; pid < 2; pid++) {
                    std::vector<int>& cur_part = partitions[pid];
                    if (cur_part.size() <= 1) continue;
                    MincutOnly::n_partitionings.fetch_add(1, std::memory_order_relaxed);
                    std::vector<std::vector<int>> ccs = GetConnectedComponentsOnPartition(
                        &induced_subgraph, cur_part);
                    for (auto& cc : ccs) {
                        MincutOnly::MinCutWorkerRecursive(
                            &induced_subgraph,
                            cc,
                            child_to_global,
                            current_connectedness_criterion,
                            connectedness_criterion_c, connectedness_criterion_x,
                            pre_computed_log,
                            connectedness_criterion_custom_string,
                            mincut_type);
                    }
                }
            } else {
                MincutOnly::n_well_connected.fetch_add(1, std::memory_order_relaxed);
                std::vector<int> done_cluster(child_to_global.begin(), child_to_global.end());
                MincutOnly::done_being_mincut_clusters.push(done_cluster);
            }

            igraph_vector_int_destroy(&nodes_to_keep);
            igraph_vector_int_destroy(&new_id_to_parent_id_vec);
            igraph_destroy(&induced_subgraph);
        }

        static inline void MinCutWorker(igraph_t* graph, ConnectednessCriterion current_connectedness_criterion, double connectedness_criterion_c, double connectedness_criterion_x, double pre_computed_log, const std::string& connectedness_criterion_custom_string, std::string mincut_type = "cactus") {
            while (true) {
                std::unique_lock<std::mutex> to_be_mincut_lock{to_be_mincut_mutex};
                to_be_mincut_condition_variable.wait(to_be_mincut_lock, []() {
                    return !MincutOnly::to_be_mincut_clusters.empty();
                });
                std::vector<int> current_cluster = MincutOnly::to_be_mincut_clusters.front();
                MincutOnly::to_be_mincut_clusters.pop();
                to_be_mincut_lock.unlock();
                if(current_cluster.size() == 1 || current_cluster[0] == -1) {
                    // done with work!
                    /* std::cerr << "thread done" << std::endl; */
                    return;
                }
                /* std::cerr << "processing cluster of size:" << std::to_string(current_cluster.size()) << std::endl; */
                /* std::cerr << "current cluster size: " << current_cluster.size() << std::endl; */
                igraph_vector_int_t nodes_to_keep;
                igraph_vector_int_t new_id_to_old_id_map;
                igraph_vector_int_init(&new_id_to_old_id_map, current_cluster.size());
                igraph_vector_int_init(&nodes_to_keep, current_cluster.size());
                for(size_t i = 0; i < current_cluster.size(); i ++) {
                    /* std::cerr << "node to keep[i]=" << std::to_string(current_cluster.at(i)) << std::endl; */
                    VECTOR(nodes_to_keep)[i] = current_cluster[i];
                }
                igraph_t induced_subgraph;
                // technically could just pass in the nodes and edges info directly by iterating through the edges and checking if it's inter vs intracluster
                // likely not too much of a memory or time overhead
                igraph_induced_subgraph_map(graph, &induced_subgraph, igraph_vss_vector(&nodes_to_keep), IGRAPH_SUBGRAPH_CREATE_FROM_SCRATCH, NULL, &new_id_to_old_id_map);
                /* std::cerr << "induced subgraph" << std::endl; */

                /* igraph_vector_int_t node_degrees; */
                /* igraph_vector_int_init(&node_degrees, 0); */
                /* igraph_degree(induced_subgraph, node_degrees, igraph_vss_vector(&nodes_to_keep), IGRAPH_ALL, IGRAPH_NO_LOOPS); */
                /* for(size_t i = 0; i < ; i ++) { */
                /*     VECTOR(nodes_to_keep)[i] = current_cluster[i]; */
                /* } */
                /* while(!ConstrainedClustering::IsWellConnected */

                /* SetIgraphAllEdgesWeight(&induced_subgraph, 1.0); */
                MinCutCustom mcc(&induced_subgraph, mincut_type);
                int edge_cut_size = mcc.ComputeMinCut();
                std::vector<int> in_partition = mcc.GetInPartition();
                std::vector<int> out_partition = mcc.GetOutPartition();
                /* std::cerr << "got the cuts into " << std::to_string(in_partition.size()) << " and " << std::to_string(out_partition.size()) << std::endl; */
                bool current_criterion = ConstrainedClustering::IsWellConnected(current_connectedness_criterion, connectedness_criterion_c, connectedness_criterion_x, pre_computed_log, in_partition.size(), out_partition.size(), edge_cut_size, connectedness_criterion_custom_string);
                /* if(connectedness_criterion == ConnectednessCriterion::Simple) { */
                /*     current_criterion = ConstrainedClustering::IsConnected(edge_cut_size); */
                /* } else if(connectedness_criterion == ConnectednessCriterion::Well) { */
                /*     current_criterion = ConstrainedClustering::IsWellConnected(in_partition, out_partition, edge_cut_size, &induced_subgraph); */
                /* } */

                if(!current_criterion) {
                    /* for(size_t i = 0; i < in_partition.size(); i ++) { */
                    /*     in_partition[i] = VECTOR(new_id_to_old_id_map)[in_partition[i]]; */
                    /* } */
                    /* for(size_t i = 0; i < out_partition.size(); i ++) { */
                    /*     out_partition[i] = VECTOR(new_id_to_old_id_map)[out_partition[i]]; */
                    /* } */
                    if(in_partition.size() > 1) {
                        std::vector<std::vector<int>> in_clusters = GetConnectedComponentsOnPartition(&induced_subgraph, in_partition);
                        for(size_t i = 0; i < in_clusters.size(); i ++) {
                            std::vector<int> translated_in_clusters;
                            for(size_t j = 0; j < in_clusters[i].size(); j ++) {
                                translated_in_clusters.push_back(VECTOR(new_id_to_old_id_map)[in_clusters[i][j]]);
                            }
                            {
                                std::lock_guard<std::mutex> to_be_mincut_guard(MincutOnly::to_be_mincut_mutex);
                                MincutOnly::to_be_mincut_clusters.push(translated_in_clusters);
                            }
                        }
                    }
                    if(out_partition.size() > 1) {
                        std::vector<std::vector<int>> out_clusters = GetConnectedComponentsOnPartition(&induced_subgraph, out_partition);
                        for(size_t i = 0; i < out_clusters.size(); i ++) {
                            std::vector<int> translated_out_clusters;
                            for(size_t j = 0; j < out_clusters[i].size(); j ++) {
                                translated_out_clusters.push_back(VECTOR(new_id_to_old_id_map)[out_clusters[i][j]]);
                            }
                            {
                                std::lock_guard<std::mutex> to_be_mincut_guard(MincutOnly::to_be_mincut_mutex);
                                MincutOnly::to_be_mincut_clusters.push(translated_out_clusters);
                            }
                        }
                    }
                } else {
                    std::unique_lock<std::mutex> done_being_mincut_lock{MincutOnly::done_being_mincut_mutex};
                    MincutOnly::done_being_mincut_clusters.push(current_cluster);
                    done_being_mincut_lock.unlock();
                }
                igraph_vector_int_destroy(&nodes_to_keep);
                igraph_vector_int_destroy(&new_id_to_old_id_map);
                igraph_destroy(&induced_subgraph);
            }
        }
    private:
        static inline std::mutex to_be_mincut_mutex;
        static inline std::condition_variable to_be_mincut_condition_variable;
        static inline std::queue<std::vector<int>> to_be_mincut_clusters;
        static inline std::mutex done_being_mincut_mutex;
        static inline std::queue<std::vector<int>> done_being_mincut_clusters;
};

#endif
