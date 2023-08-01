#pragma once
#ifndef LION_H
#define LION_H

#include "shortest_path.h"
#include "KAF.h"

void init_LION(model& our_model);
void lixel_augmentation(model& our_model);
void lixel_aggregation(model& our_model);
void init_one_D_KDV(model& our_model);
void one_D_KDV(model& our_model);

#endif