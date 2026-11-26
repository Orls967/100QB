#include "common.h"

Config load_config(const char *path) {
    Config c; memset(&c, 0, sizeof(c));
    FILE *f = fopen(path, "r");
    if (!f) { fprintf(stderr, "cannot open %s\n", path); exit(1); }
    char line[256];
    while (fgets(line, sizeof(line), f)) {
        char *eq = strchr(line, '=');
        if (!eq) continue;
        *eq = '\0';
        char *k = line, *v = eq + 1;
        while (*k == ' ' || *k == '\t') k++;
        int64_t val = strtoll(v, NULL, 10);
        if (!strcmp(k, "master_seed")) c.master_seed = val;
        else if (!strcmp(k, "slots_per_day")) c.slots_per_day = val;
        else if (!strcmp(k, "slots_per_hour")) c.slots_per_hour = val;
        else if (!strcmp(k, "parity_days")) c.parity_days = val;
        else if (!strcmp(k, "arr_scale_bp")) c.arr_scale_bp = val;
        else if (!strcmp(k, "checkin_servers")) c.srv_cnt[STN_CHECKIN] = val;
        else if (!strcmp(k, "security_servers")) c.srv_cnt[STN_SECURITY] = val;
        else if (!strcmp(k, "imigrasi_servers")) c.srv_cnt[STN_IMIGRASI] = val;
        else if (!strcmp(k, "gate_servers")) c.srv_cnt[STN_GATE] = val;
        else if (!strcmp(k, "belt_servers")) c.srv_cnt[STN_BELT] = val;
        else if (!strcmp(k, "customs_servers")) c.srv_cnt[STN_CUSTOMS] = val;
        else if (!strcmp(k, "checkin_mu_bp")) c.mu_bp[STN_CHECKIN] = val;
        else if (!strcmp(k, "security_mu_bp")) c.mu_bp[STN_SECURITY] = val;
        else if (!strcmp(k, "security_extra_prob_bp")) c.security_extra_prob_bp = val;
        else if (!strcmp(k, "security_extra_mu_bp")) c.security_extra_mu_bp = val;
        else if (!strcmp(k, "imigrasi_mu_bp")) c.mu_bp[STN_IMIGRASI] = val;
        else if (!strcmp(k, "gate_mu_bp")) c.mu_bp[STN_GATE] = val;
        else if (!strcmp(k, "belt_mu_bp")) c.mu_bp[STN_BELT] = val;
        else if (!strcmp(k, "customs_mu_bp")) c.mu_bp[STN_CUSTOMS] = val;
        else if (!strcmp(k, "belt_cap")) c.belt_cap = val;
        else if (!strcmp(k, "break_prob_bp")) c.break_prob_bp = val;
        else if (!strcmp(k, "repair_prob_bp")) c.repair_prob_bp = val;
        else if (!strcmp(k, "patience_slots")) c.patience_slots = val;
        else if (!strcmp(k, "gate_close_lead")) c.gate_close_lead = val;
        else if (!strcmp(k, "bag_miss_bp")) c.bag_miss_bp = val;
        else if (!strcmp(k, "intl_prob_bp")) c.intl_prob_bp = val;
        else if (!strcmp(k, "biz_prob_bp")) c.biz_prob_bp = val;
        else if (!strcmp(k, "transit_prob_bp")) c.transit_prob_bp = val;
        else if (!strcmp(k, "rebook_prob_bp")) c.rebook_prob_bp = val;
    }
    fclose(f);
    return c;
}

RunParams default_rp(const Config *c, int64_t days) {
    RunParams rp; memset(&rp, 0, sizeof(rp));
    rp.days = days;
    rp.arr_scale_bp = c->arr_scale_bp;
    rp.break_prob_bp = c->break_prob_bp;
    rp.repair_prob_bp = c->repair_prob_bp;
    rp.belt_cap = c->belt_cap;
    rp.gate_close_lead = c->gate_close_lead;
    for (int i = 0; i < N_STN; ++i) { rp.srv_ovr[i] = -1; rp.mu_ovr[i] = -1; }
    rp.priority_on = 1; rp.breakdown_on = 1; rp.roster_on = 1; rp.customs_on = 0;
    rp.audit_on = 0;
    return rp;
}
