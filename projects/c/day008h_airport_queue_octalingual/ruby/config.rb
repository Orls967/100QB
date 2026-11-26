# frozen_string_literal: true

SLOTS_PER_DAY = 8640
SLOTS_PER_HOUR = 360
N_STN = 6
MAX_SRV = 16
MAX_Q = 4096
MAX_POOL = 8192
AUDIT_DAYS = 3

STN_CHECKIN = 0; STN_SECURITY = 1; STN_IMIGRASI = 2; STN_GATE = 3; STN_BELT = 4; STN_CUSTOMS = 5
CLS_ECON = 0; CLS_BIZ = 1; CLS_TRANSIT = 2
DEST_DOM = 0; DEST_INTL = 1
SRV_IDLE = 0; SRV_BUSY = 1; SRV_DOWN = 2

RunParams = Struct.new(
  :days, :arr_scale_bp, :break_prob_bp, :repair_prob_bp, :belt_cap, :gate_close_lead,
  :srv_ovr, :mu_ovr, :priority_on, :breakdown_on, :roster_on, :customs_on,
  :permute_srv, :timescale_div2, :audit_on,
  keyword_init: true
)

def default_rp(c, days)
  RunParams.new(
    days: days, arr_scale_bp: c["arr_scale_bp"], break_prob_bp: c["break_prob_bp"],
    repair_prob_bp: c["repair_prob_bp"], belt_cap: c["belt_cap"], gate_close_lead: c["gate_close_lead"],
    srv_ovr: [-1] * N_STN, mu_ovr: [-1] * N_STN,
    priority_on: 1, breakdown_on: 1, roster_on: 1, customs_on: 0, permute_srv: 0, timescale_div2: 0, audit_on: 0
  )
end

def load_config(path)
  c = {}
  File.readlines(path, chomp: true).each do |line|
    next unless line.include?("=")

    k, v = line.split("=", 2)
    c[k.strip] = v.strip.to_i
  end
  c["srv_cnt"] = [c["checkin_servers"], c["security_servers"], c["imigrasi_servers"], c["gate_servers"], c["belt_servers"], c["customs_servers"]]
  c["mu_bp"] = [c["checkin_mu_bp"], c["security_mu_bp"], c["imigrasi_mu_bp"], c["gate_mu_bp"], c["belt_mu_bp"], c["customs_mu_bp"]]
  c
end
