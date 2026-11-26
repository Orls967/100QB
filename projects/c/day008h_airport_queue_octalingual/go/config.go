package main

import (
	"os"
	"strconv"
	"strings"
)

const (
	SlotsPerDay  int64 = 8640
	SlotsPerHour int64 = 360
	NStn               = 6
	MaxSrv             = 16
	MaxQ               = 4096
	MaxPool            = 8192
	AuditDays    int64 = 3

	StnCheckin  = 0
	StnSecurity = 1
	StnImigrasi = 2
	StnGate     = 3
	StnBelt     = 4
	StnCustoms  = 5

	ClsEcon    int32 = 0
	ClsBiz     int32 = 1
	ClsTransit int32 = 2
	DestDom    int32 = 0
	DestIntl   int32 = 1
	SrvIdle    int32 = 0
	SrvBusy    int32 = 1
	SrvDown    int32 = 2
)

type Config struct {
	masterSeed, slotsPerDay, slotsPerHour, parityDays, arrScaleBp int64
	srvCnt                                                        [NStn]int64
	muBp                                                          [NStn]int64
	securityExtraProbBp, securityExtraMuBp                        int64
	beltCap, breakProbBp, repairProbBp, patienceSlots             int64
	gateCloseLead, bagMissBp, intlProbBp, bizProbBp               int64
	transitProbBp, rebookProbBp                                   int64
}

type RunParams struct {
	days, arrScaleBp, breakProbBp, repairProbBp, beltCap, gateCloseLead int64
	srvOvr                                                              [NStn]int64
	muOvr                                                               [NStn]int64
	priorityOn, breakdownOn, rosterOn, customsOn                        int32
	permuteSrv, timescaleDiv2, auditOn                                  int32
}

func defaultRP(c Config, days int64) RunParams {
	return RunParams{
		days: days, arrScaleBp: c.arrScaleBp, breakProbBp: c.breakProbBp, repairProbBp: c.repairProbBp,
		beltCap: c.beltCap, gateCloseLead: c.gateCloseLead,
		srvOvr: [NStn]int64{-1, -1, -1, -1, -1, -1}, muOvr: [NStn]int64{-1, -1, -1, -1, -1, -1},
		priorityOn: 1, breakdownOn: 1, rosterOn: 1,
	}
}

func loadConfig(path string) Config {
	data, err := os.ReadFile(path)
	if err != nil {
		panic(err)
	}
	var c Config
	for _, line := range strings.Split(string(data), "\n") {
		k, v, ok := strings.Cut(line, "=")
		if !ok {
			continue
		}
		val, _ := strconv.ParseInt(strings.TrimSpace(v), 10, 64)
		switch strings.TrimSpace(k) {
		case "master_seed":
			c.masterSeed = val
		case "slots_per_day":
			c.slotsPerDay = val
		case "slots_per_hour":
			c.slotsPerHour = val
		case "parity_days":
			c.parityDays = val
		case "arr_scale_bp":
			c.arrScaleBp = val
		case "checkin_servers":
			c.srvCnt[0] = val
		case "security_servers":
			c.srvCnt[1] = val
		case "imigrasi_servers":
			c.srvCnt[2] = val
		case "gate_servers":
			c.srvCnt[3] = val
		case "belt_servers":
			c.srvCnt[4] = val
		case "customs_servers":
			c.srvCnt[5] = val
		case "checkin_mu_bp":
			c.muBp[0] = val
		case "security_mu_bp":
			c.muBp[1] = val
		case "security_extra_prob_bp":
			c.securityExtraProbBp = val
		case "security_extra_mu_bp":
			c.securityExtraMuBp = val
		case "imigrasi_mu_bp":
			c.muBp[2] = val
		case "gate_mu_bp":
			c.muBp[3] = val
		case "belt_mu_bp":
			c.muBp[4] = val
		case "customs_mu_bp":
			c.muBp[5] = val
		case "belt_cap":
			c.beltCap = val
		case "break_prob_bp":
			c.breakProbBp = val
		case "repair_prob_bp":
			c.repairProbBp = val
		case "patience_slots":
			c.patienceSlots = val
		case "gate_close_lead":
			c.gateCloseLead = val
		case "bag_miss_bp":
			c.bagMissBp = val
		case "intl_prob_bp":
			c.intlProbBp = val
		case "biz_prob_bp":
			c.bizProbBp = val
		case "transit_prob_bp":
			c.transitProbBp = val
		case "rebook_prob_bp":
			c.rebookProbBp = val
		}
	}
	return c
}
