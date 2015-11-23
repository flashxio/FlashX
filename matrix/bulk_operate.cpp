/*
 * Copyright 2015 Open Connectome Project (http://openconnecto.me)
 * Written by Da Zheng (zhengda1936@gmail.com)
 *
 * This file is part of FlashMatrix.
 *
 * Licensed under the Apache License, Version 2.0 (the "License");
 * you may not use this file except in compliance with the License.
 * You may obtain a copy of the License at
 *
 *     http://www.apache.org/licenses/LICENSE-2.0
 *
 * Unless required by applicable law or agreed to in writing, software
 * distributed under the License is distributed on an "AS IS" BASIS,
 * WITHOUT WARRANTIES OR CONDITIONS OF ANY KIND, either express or implied.
 * See the License for the specific language governing permissions and
 * limitations under the License.
 */
#include "log.h"

#include "bulk_operate.h"
#include "bulk_operate_ext.h"

namespace fm
{

agg_operate::const_ptr agg_operate::create(bulk_operate::const_ptr agg,
			bulk_operate::const_ptr combine)
{
	if (agg == NULL) {
		BOOST_LOG_TRIVIAL(error)
			<< "the agg operator is required in agg_operate";
		return const_ptr();
	}
	if (combine && agg->get_output_type() != combine->get_left_type()
			&& combine->get_left_type() != combine->get_output_type()) {
		BOOST_LOG_TRIVIAL(error)
			<< "the agg and combine operators have incompatible types";
		return const_ptr();
	}
	return const_ptr(new agg_operate(agg, combine));
}

agg_operate::const_ptr agg_operate::create(bulk_operate::const_ptr agg)
{
	if (agg->get_left_type() != agg->get_right_type()
			|| agg->get_left_type() != agg->get_output_type()) {
		BOOST_LOG_TRIVIAL(error)
			<< "The agg operator needs to have the same input and output type";
		return const_ptr();
	}
	return const_ptr(new agg_operate(agg, agg));
}

}
