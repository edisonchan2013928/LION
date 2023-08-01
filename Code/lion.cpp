#include "lion.h"

void init_LION(model& our_model)
{
	int num_deg;

	if (our_model.k_type == 1) //Triangular kernel
		our_model.num_deg = 2;
	if (our_model.k_type == 2) //Epanechnikov kernel
		our_model.num_deg = 3;
	if (our_model.k_type == 3) //Quartic kernel
		our_model.num_deg = 5;

	num_deg = our_model.num_deg;
	for (int e = 0; e < our_model.m; e++)
	{
		Edge& edge = our_model.edge_set[e];

		for (int q = 0; q < (int)edge.lixel_vec.size(); q++)
		{
			//degree = 3 (0, 1, and 2)
			edge.lixel_vec[q].alpha_B_x_deg_vec = new double[num_deg];
			edge.lixel_vec[q].alpha_B_y_deg_vec = new double[num_deg];
			edge.lixel_vec[q].alpha_R_x_deg_vec = new double[num_deg];
			edge.lixel_vec[q].alpha_R_y_deg_vec = new double[num_deg];
		}
	}

	our_model.ell.a_L = new double[our_model.num_deg];
	our_model.ell.a_R = new double[our_model.num_deg];
}

void lixel_augmentation(model& our_model)
{
	int num_deg = our_model.num_deg;
	int e_index;
	int node_index_pi_a; int node_index_pi_b;
	int node_index_pi_u; int node_index_pi_v;
	double temp_value;
	double dist_a_p;
	double dist_b_p;
	int q_r;
	int q_l;

	our_model.lixel_total_number = 0;
	for (int wide_e_index = 0; wide_e_index < our_model.m; wide_e_index++)
	{
		Edge& wide_edge = our_model.edge_set[wide_e_index];
		our_model.lixel_total_number += (int)wide_edge.lixel_vec.size();

		for (int q = 0; q < (int)wide_edge.lixel_vec.size(); q++)
		{
			for (int d = 0; d < num_deg; d++)
			{
				wide_edge.lixel_vec[q].alpha_B_x_deg_vec[d] = 0;
				wide_edge.lixel_vec[q].alpha_B_y_deg_vec[d] = 0;
				wide_edge.lixel_vec[q].alpha_R_x_deg_vec[d] = 0;
				wide_edge.lixel_vec[q].alpha_R_y_deg_vec[d] = 0;
			}
		}

		our_model.cur_edge_index = wide_e_index;

		//Obtain SPD from node x (node a) to all other nodes.
		clear_scan_edge_list(our_model);
		our_model.sel_node_index = wide_edge.n1;
		our_model.node_index_pi_a = our_model.sel_node_index;
		dijkstra(our_model);
		copy_sp_info(our_model, true);

		//Obtain SPD from node y (node b) to all other nodes.
		our_model.sel_node_index = wide_edge.n2;
		our_model.node_index_pi_b = our_model.sel_node_index;
		dijkstra(our_model);
		copy_sp_info(our_model, false);

		node_index_pi_a = our_model.node_index_pi_a;
		node_index_pi_b = our_model.node_index_pi_b;
		for (int e = 0; e < (int)our_model.access_edge_list.size(); e++)
		{
			e_index = our_model.access_edge_list[e];

			if (e_index == wide_e_index)
				continue;
			
			Edge& edge = our_model.edge_set[e_index];
			node_index_pi_u = edge.n1;
			node_index_pi_v = edge.n2;

			for (int p = 0; p < (int)edge.PS.size(); p++)
			{
				dist_a_p = min(our_model.sp_node_vec_node_a[node_index_pi_u].cur_sp_value + edge.PS[p].dist_n1,
					our_model.sp_node_vec_node_a[node_index_pi_v].cur_sp_value + edge.PS[p].dist_n2);
				dist_b_p = min(our_model.sp_node_vec_node_b[node_index_pi_u].cur_sp_value + edge.PS[p].dist_n1,
					our_model.sp_node_vec_node_b[node_index_pi_v].cur_sp_value + edge.PS[p].dist_n2);

				if (dist_a_p > our_model.bandwidth && dist_b_p > our_model.bandwidth) //Case 1
					continue;

				if (dist_a_p <= our_model.bandwidth && dist_b_p > our_model.bandwidth) //Case 2
				{
					//q_r = (int)floor((our_model.bandwidth - dist_a_p - wide_edge.lixel_vec[0].dist_n1) / our_model.lixel_reg_length) - 1;
					q_r = (int)floor((our_model.bandwidth - dist_a_p - wide_edge.lixel_vec[0].dist_n1) / our_model.lixel_reg_length);

					if (q_r <= -1)
						continue;

					temp_value = 1;
					for (int d = 0; d < our_model.num_deg; d++)
					{
						wide_edge.lixel_vec[q_r].alpha_B_x_deg_vec[d] += temp_value;
						temp_value *= dist_a_p;
					}
				}

				if (dist_a_p > our_model.bandwidth && dist_b_p <= our_model.bandwidth) //Case 3
				{
					q_l = wide_edge.lixel_vec.size() - (int)floor((our_model.bandwidth - dist_b_p - wide_edge.lixel_vec[wide_edge.lixel_vec.size() - 1].dist_n2) / our_model.lixel_reg_length) - 1;

					if (q_l >= wide_edge.lixel_vec.size())
						continue;

					temp_value = 1;
					for (int d = 0; d < our_model.num_deg; d++)
					{
						wide_edge.lixel_vec[q_l].alpha_B_y_deg_vec[d] += temp_value;
						temp_value *= dist_b_p;
					}
				}

				if (dist_a_p <= our_model.bandwidth && dist_b_p <= our_model.bandwidth) //Case 4
				{
					if (2 * our_model.bandwidth - dist_a_p - dist_b_p >= wide_edge.length)
					{
						//q_r = (int)floor(((wide_edge.length - dist_a_p + dist_b_p) / 2.0 - wide_edge.lixel_vec[0].dist_n1) / our_model.lixel_reg_length) - 1;
						q_r = (int)floor(((wide_edge.length - dist_a_p + dist_b_p) / 2.0 - wide_edge.lixel_vec[0].dist_n1) / our_model.lixel_reg_length);
						q_l = q_r + 1;

						if (q_r > -1 && q_l < wide_edge.lixel_vec.size())
						{
							temp_value = 1;
							for (int d = 0; d < our_model.num_deg; d++)
							{
								wide_edge.lixel_vec[q_r].alpha_B_x_deg_vec[d] += temp_value;
								temp_value *= dist_a_p;
							}

							temp_value = 1;
							for (int d = 0; d < our_model.num_deg; d++)
							{
								wide_edge.lixel_vec[q_l].alpha_B_y_deg_vec[d] += temp_value;
								temp_value *= dist_b_p;
							}
						}
						else
						{
							if (q_r <= -1)
							{
								temp_value = 1;
								for (int d = 0; d < our_model.num_deg; d++)
								{
									wide_edge.lixel_vec[q_l].alpha_B_y_deg_vec[d] += temp_value;
									temp_value *= dist_b_p;
								}
							}
							
							if (q_l >= wide_edge.lixel_vec.size())
							{
								temp_value = 1;
								for (int d = 0; d < our_model.num_deg; d++)
								{
									wide_edge.lixel_vec[q_r].alpha_B_x_deg_vec[d] += temp_value;
									temp_value *= dist_a_p;
								}
							}
						}
					}
					else
					{
						//q_r = (int)floor((our_model.bandwidth - dist_a_p - wide_edge.lixel_vec[0].dist_n1) / our_model.lixel_reg_length) - 1;
						q_r = (int)floor((our_model.bandwidth - dist_a_p - wide_edge.lixel_vec[0].dist_n1) / our_model.lixel_reg_length);

						if (q_r > -1)
						{
							temp_value = 1;
							for (int d = 0; d < our_model.num_deg; d++)
							{
								wide_edge.lixel_vec[q_r].alpha_B_x_deg_vec[d] += temp_value;
								temp_value *= dist_a_p;
							}
						}

						//q_l = wide_edge.lixel_vec.size() - (int)floor((our_model.bandwidth - dist_b_p - wide_edge.lixel_vec[wide_edge.lixel_vec.size() - 1].dist_n2) / our_model.lixel_reg_length);
						q_l = wide_edge.lixel_vec.size() - (int)floor((our_model.bandwidth - dist_b_p - wide_edge.lixel_vec[wide_edge.lixel_vec.size() - 1].dist_n2) / our_model.lixel_reg_length) - 1;

						if (q_l < wide_edge.lixel_vec.size())
						{
							temp_value = 1;
							for (int d = 0; d < our_model.num_deg; d++)
							{
								wide_edge.lixel_vec[q_l].alpha_B_y_deg_vec[d] += temp_value;
								temp_value *= dist_b_p;
							}
						}
					}
				}

			}
		}
	}
}

void lixel_aggregation(model& our_model)
{
	double*alpha_left_vec = new double[our_model.num_deg];
	double*alpha_right_vec = new double[our_model.num_deg];

	for (int wide_e_index = 0; wide_e_index < our_model.m; wide_e_index++)
	{
		Edge& wide_edge = our_model.edge_set[wide_e_index];
		for (int d = 0; d < our_model.num_deg; d++)
		{
			alpha_left_vec[d] = 0;
			alpha_right_vec[d] = 0;
		}

		for (int q = wide_edge.lixel_vec.size() - 1; q >= 0; q--)
		{
			for (int d = 0; d < our_model.num_deg; d++)
			{
				alpha_left_vec[d] += wide_edge.lixel_vec[q].alpha_B_x_deg_vec[d];
				wide_edge.lixel_vec[q].alpha_R_x_deg_vec[d] = alpha_left_vec[d];
			}
		}

		for (int q = 0; q < (int)wide_edge.lixel_vec.size(); q++)
		{
			for (int d = 0; d < our_model.num_deg; d++)
			{
				alpha_right_vec[d] += wide_edge.lixel_vec[q].alpha_B_y_deg_vec[d];
				wide_edge.lixel_vec[q].alpha_R_y_deg_vec[d] = alpha_right_vec[d];
			}
		}

		for (int q = 0; q < (int)wide_edge.lixel_vec.size(); q++)
		{
			if (our_model.k_type == 2) //Epanechnikov kernel
			{
				wide_edge.lixel_vec[q].KDE_value = (1 - (wide_edge.lixel_vec[q].dist_n1 * wide_edge.lixel_vec[q].dist_n1) / (our_model.bandwidth * our_model.bandwidth)) * wide_edge.lixel_vec[q].alpha_R_x_deg_vec[0]
					- (2.0 * wide_edge.lixel_vec[q].dist_n1 * wide_edge.lixel_vec[q].alpha_R_x_deg_vec[1]) / (our_model.bandwidth*our_model.bandwidth)
					- (wide_edge.lixel_vec[q].alpha_R_x_deg_vec[2] / (our_model.bandwidth * our_model.bandwidth))
					+ (1 - (wide_edge.lixel_vec[q].dist_n2 * wide_edge.lixel_vec[q].dist_n2) / (our_model.bandwidth * our_model.bandwidth)) * wide_edge.lixel_vec[q].alpha_R_y_deg_vec[0]
					- (2.0 * wide_edge.lixel_vec[q].dist_n2 * wide_edge.lixel_vec[q].alpha_R_y_deg_vec[1]) / (our_model.bandwidth*our_model.bandwidth)
					- (wide_edge.lixel_vec[q].alpha_R_y_deg_vec[2] / (our_model.bandwidth * our_model.bandwidth));
			}

			//Other kernels (Code here)
		}
	}

	delete[] alpha_left_vec;
	delete[] alpha_right_vec;
}

void init_one_D_KDV(model& our_model)
{
	for (int d = 0; d < our_model.num_deg; d++)
	{
		our_model.ell.a_L[d] = 0;
		our_model.ell.a_R[d] = 0;
	}
}

void obtain_sorted_end_vec(model& our_model, Edge& edge, vector<end_point>& sorted_end_vec)
{
	vector<end_point> end_vec_left;
	vector<end_point> end_vec_right;
	end_point temp_e_point;
	int end_point_num = 2 * edge.PS.size();
	int left_num = 0;
	int right_num = 0;

	for (int p = 0; p < (int)edge.PS.size(); p++)
	{
		temp_e_point.index = p;
		temp_e_point.isLeft = true;
		temp_e_point.value = edge.PS[p].dist_n1 - our_model.bandwidth;
		
		end_vec_left.push_back(temp_e_point);

		temp_e_point.isLeft = false;
		temp_e_point.value = edge.PS[p].dist_n1 + our_model.bandwidth;
		
		end_vec_right.push_back(temp_e_point);
	}

	for (int i = 0; i < end_point_num; i++)
	{
		if (left_num == (int)edge.PS.size())
		{
			sorted_end_vec.push_back(end_vec_right[right_num]);
			right_num++;
			continue;
		}
		
		if (right_num == (int)edge.PS.size())
		{
			sorted_end_vec.push_back(end_vec_left[left_num]);
			left_num++;
			continue;
		}

		if (end_vec_left[left_num].value < end_vec_right[right_num].value)
		{
			sorted_end_vec.push_back(end_vec_left[left_num]);
			left_num++;
		}
		else
		{
			sorted_end_vec.push_back(end_vec_right[right_num]);
			right_num++;
		}
	}
}

void sweep_line_algorithm(model& our_model, Edge& edge, vector<end_point>& sorted_end_vec)
{
	int lixel_size = edge.lixel_vec.size();
	int end_point_num = 2 * edge.PS.size();
	int lixel_index = 0;
	int end_point_index = 0;
	double temp_value;

	for (int i = 0; i < lixel_size + end_point_num; i++)
	{
		if (lixel_index == lixel_size || end_point_index == end_point_num)
			break;

		//if (edge.lixel_vec[lixel_index].dist_n1 < edge.PS[sorted_end_vec[end_point_index].index].dist_n1)
		if (edge.lixel_vec[lixel_index].dist_n1 < sorted_end_vec[end_point_index].value)
		{
			if (our_model.k_type == 2) //Epanechnikov kernel
			{
				edge.lixel_vec[lixel_index].KDE_value += (1.0 - (edge.lixel_vec[lixel_index].dist_n1 * edge.lixel_vec[lixel_index].dist_n1) / (our_model.bandwidth * our_model.bandwidth)) * our_model.ell.a_L[0]
					+ (2.0 * edge.lixel_vec[lixel_index].dist_n1 * our_model.ell.a_L[1]) / (our_model.bandwidth * our_model.bandwidth)
					- our_model.ell.a_L[2] / (our_model.bandwidth * our_model.bandwidth)
					- (1.0 - (edge.lixel_vec[lixel_index].dist_n1 * edge.lixel_vec[lixel_index].dist_n1) / (our_model.bandwidth * our_model.bandwidth)) * our_model.ell.a_R[0]
					- (2.0 * edge.lixel_vec[lixel_index].dist_n1 * our_model.ell.a_R[1]) / (our_model.bandwidth * our_model.bandwidth)
					+ our_model.ell.a_R[2] / (our_model.bandwidth * our_model.bandwidth);
			}

			//Other kernels (Code here)

			lixel_index++;
		}
		else
		{
			if (sorted_end_vec[end_point_index].isLeft == true)
			{
				temp_value = 1;
				for (int d = 0; d < our_model.num_deg; d++)
				{
					our_model.ell.a_L[d] += temp_value;
					temp_value *= edge.PS[sorted_end_vec[end_point_index].index].dist_n1;
				}
			}
			else
			{
				temp_value = 1;
				for (int d = 0; d < our_model.num_deg; d++)
				{
					our_model.ell.a_R[d] += temp_value;
					temp_value *= edge.PS[sorted_end_vec[end_point_index].index].dist_n1;
				}
			}

			end_point_index++;
		}
	}
}

void one_D_KDV(model& our_model) //Section 3.1
{
	vector<end_point> sorted_end_vec;

	for (int wide_e = 0; wide_e < our_model.m; wide_e++)
	{
		init_one_D_KDV(our_model);
		Edge& edge = our_model.edge_set[wide_e];

		if (edge.PS.size() == 0)
			continue;

		obtain_sorted_end_vec(our_model, edge, sorted_end_vec);
		sweep_line_algorithm(our_model, edge, sorted_end_vec);
		sorted_end_vec.clear();
	}
}