#ifndef __FG2TPETRA_H__
#define __FG2TPETRA_H__

#include "sparse_matrix.h"
#include "vertex.h"
#include "vertex_index.h"
#include "in_mem_storage.h"
#include "io_interface.h"

#include "Tpetra_Map.hpp"
#include "Tpetra_MultiVector.hpp"
#include "Tpetra_DefaultPlatform.hpp"
#include "Tpetra_CrsMatrix_decl.hpp"
#include "Tpetra_CrsMatrix_def.hpp"

using Teuchos::RCP;
using Teuchos::rcp;
using Teuchos::ArrayRCP;
using std::cerr;
using std::cout;
using std::endl;

typedef size_t local_ordinal_type;
typedef size_t global_ordinal_type;

typedef Tpetra::MultiVector<double, global_ordinal_type, global_ordinal_type> MV;

typedef Tpetra::Map<local_ordinal_type, global_ordinal_type> map_type;
typedef Tpetra::CrsMatrix<double, local_ordinal_type, global_ordinal_type> crs_matrix_type;

RCP<crs_matrix_type> create_crs(const std::string &graph_file,
		fg::vertex_index::ptr index, fg::edge_type type,
		RCP<map_type> map, int my_rank, size_t edge_data_size);

#endif
