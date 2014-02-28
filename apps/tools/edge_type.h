#ifndef __EDGE_TYPE_H__
#define __EDGE_TYPE_H__

/**
 * The type of edge data.
 */
enum {
	DEFAULT_TYPE,
	EDGE_COUNT,
	EDGE_TIMESTAMP,
};

struct str2int_pair {
	std::string str;
	int number;
};
static struct str2int_pair edge_type_map[] = {
	{"count", EDGE_COUNT},
	{"timestamp", EDGE_TIMESTAMP},
};
static int type_map_size = sizeof(edge_type_map) / sizeof(edge_type_map[0]);

static inline int conv_edge_type_str2int(const std::string &type_str)
{
	for (int i = 0; i < type_map_size; i++) {
		if (edge_type_map[i].str == type_str) {
			return edge_type_map[i].number;
		}
	}
	return DEFAULT_TYPE;
}

class ts_edge_data
{
	time_t timestamp;
	float weight;
public:
	ts_edge_data() {
		timestamp = 0;
		weight = 1;
	}

	ts_edge_data(time_t timestamp, float weight) {
		this->timestamp = timestamp;
		this->weight = weight;
	}

	time_t get_timestamp() const {
		return timestamp;
	}

	float get_weight() const {
		return weight;
	}
};

#endif
