
#include "pop.hpp"

void population::eval_fit(agent* a, int i) {
    
    if( g_debug )
        cout << "eval_fit(agent* a, int i) i:" << i << endl;

    a->fitness[i] = train_e.eval_fit(a->member[i]);
}

void population::recombine(agent* a, agent* b) {

    if( g_debug )
        cout << "recombine(agent* a, agent* b) a->depth:" << a->depth << " a->degree" << a->degree << " b->depth" << b->depth << " b->degree" << b->degree << endl;

    switch (g_rint(1, 3)) {
        case 1:
            variable_intersect(a, b);
            break;
        case 2:
            variable_union(a, b);
            break;
        case 3:
            variable_symdiff(a, b);
        default:
            break;
    }

    eval_fit(a, 1);	//1=current
    //    //Pablo: Asked me to Run LS here instead of fitness eval in previous line
    //    optimize opt(train_data, a->member[1], objec, delta);
    //    a->fitness[1] = opt.run();

    a->update_pocket();
    a->propagate_pocket();
}

void population::mutate(agent* a) {
    
    if( g_debug )
        cout << "mutate(agent* a) a->depth:" << a->depth << " a->degree:" << a->degree << endl;

    //moh: if(rint(1, 4) == 1) {
    if (g_rreal() < g_mu_rate) {
        
        // a rather major mutation
        // only mutate if current is too close or much worse
        // that pocket
        // If depth is locked we must not hard mutate s   
        if ( (a->fitness[1] < 1.2 * a->fitness[0] || a->fitness[1] > 2 * a->fitness[0]) && g_depth_lock == -1 )
            feature_toggle(a->member[1]);
        else
            // or we perturb only slightly
            feature_mod(a->member[1]);
            
    }

    eval_fit(a, 1);
    // Pablo: Asked me to Run LS here instead of fitness eval in previous line
    //        cout<< "Start LS after Mutate" <<endl;

    //    optimize opt(train_data, a->member[1], objec, delta);
    //    a->fitness[1] = opt.run();
    // 10-Aug-2020: Mohammad removed running LS 

    a->update_pocket();
    a->propagate_pocket();
}

void population::local_search(agent* a) {
    
    if( g_debug )
        cout << "local_search(agent* a) a->depth:" << a->depth << " a->degree" << a->degree << endl;
    
    if (a->depth > 1) {
        for (size_t i = 0; i < a->degree; ++i)
            local_search(a->children[i]);
    }

    optimize opt(train_data, a->member[1], objec, delta);
    a->fitness[1] = opt.run();

    if (a->fitness[0] > a->fitness[1]) {
        optimize more_opt(train_data, a->member[0], objec, delta);
        a->fitness[0] = more_opt.run();
    }

    a->movedown_pocket();
}

