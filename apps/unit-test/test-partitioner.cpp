#include <stdlib.h>

#include "partitioner.h"

const int num_parts = 16;
const int M = 1024 * 1024;

void test_partitioner(graph_partitioner &partitioner)
{
	for (int k = 0; k < 100; k++) {
		std::vector<vertex_id_t> parts[num_parts];
		size_t num_vertices = random() % M + M;
		printf("there are %ld vertices\n", num_vertices);
		for (int i = 0; i < num_parts; i++) {
			partitioner.get_all_vertices_in_part(i, num_vertices, parts[i]);
			size_t computed_part_size = partitioner.get_part_size(i,
					num_vertices);
			assert(computed_part_size == parts[i].size());
		}
		for (vertex_id_t id = 0; id < num_vertices; id++) {
			int part_id;
			off_t off;
			partitioner.map2loc(id, part_id, off);
			assert(part_id == partitioner.map(id));
			assert(parts[part_id][off] == id);
		}
		for (int part_id = 0; part_id < num_parts; part_id++) {
			for (off_t off = 0; off < parts[part_id].size(); off++) {
				vertex_id_t id;
				partitioner.loc2map(part_id, off, id);
				assert(id == parts[part_id][off]);
			}
		}
		size_t tot = 0;
		for (int i = 0; i < num_parts; i++)
			tot += parts[i].size();
		printf("There are %ld vertices in all partitions\n", tot);
		assert(num_vertices == tot);
	}
}

int main()
{
	printf("test range_graph_partitioner\n");
	range_graph_partitioner r_partitioner(num_parts);
	test_partitioner(r_partitioner);

	printf("test modulo_graph_partitioner\n");
	modulo_graph_partitioner m_partitioner(num_parts);
	test_partitioner(m_partitioner);

}
