//! Battery storage specification: 4 MWh capacity, 0.1 MWh SoC grid (41 states),
//! eta_c = eta_d = 0.90, 5 discrete actions, degradation 2 EUR/MWh grid-side, and terminal value.

pub const NUM_SOC_STATES: usize = 41;
pub const NUM_ACTIONS: usize = 5;

pub const ACTION_CHARGE_2MW: usize = 0;
pub const ACTION_CHARGE_1MW: usize = 1;
pub const ACTION_IDLE: usize = 2;
pub const ACTION_DISCHARGE_09MW: usize = 3;
pub const ACTION_DISCHARGE_18MW: usize = 4;

#[derive(Clone, Copy, Debug)]
pub struct BatteryAction {
    pub id: usize,
    pub name: &'static str,
    pub power_mw: f64,
    pub grid_through_mw: f64,
    pub stored_delta: f64,
    pub delta_soc_idx: isize,
}

pub const ACTIONS: [BatteryAction; NUM_ACTIONS] = [
    BatteryAction {
        id: ACTION_CHARGE_2MW,
        name: "charge_2mw",
        power_mw: 2.0,
        grid_through_mw: 2.0,
        stored_delta: 1.8,
        delta_soc_idx: 18,
    },
    BatteryAction {
        id: ACTION_CHARGE_1MW,
        name: "charge_1mw",
        power_mw: 1.0,
        grid_through_mw: 1.0,
        stored_delta: 0.9,
        delta_soc_idx: 9,
    },
    BatteryAction {
        id: ACTION_IDLE,
        name: "idle",
        power_mw: 0.0,
        grid_through_mw: 0.0,
        stored_delta: 0.0,
        delta_soc_idx: 0,
    },
    BatteryAction {
        id: ACTION_DISCHARGE_09MW,
        name: "discharge_09mw",
        power_mw: -0.9,
        grid_through_mw: 0.9,
        stored_delta: -1.0,
        delta_soc_idx: -10,
    },
    BatteryAction {
        id: ACTION_DISCHARGE_18MW,
        name: "discharge_18mw",
        power_mw: -1.8,
        grid_through_mw: 1.8,
        stored_delta: -2.0,
        delta_soc_idx: -20,
    },
];

#[derive(Clone, Copy, Debug)]
pub struct BatterySpec {
    pub capacity_mwh: f64,
    pub initial_soc_mwh: f64,
    pub initial_soc_idx: usize,
    pub soc_step_mwh: f64,
    pub num_soc_states: usize,
    pub eta_c: f64,
    pub eta_d: f64,
    pub max_charge_mw: f64,
    pub max_discharge_mw: f64,
    pub degradation_cost: f64,
    pub terminal_ref_price: f64,
}

impl BatterySpec {
    pub fn default_spec() -> Self {
        Self {
            capacity_mwh: 4.0,
            initial_soc_mwh: 2.0,
            initial_soc_idx: 20,
            soc_step_mwh: 0.1,
            num_soc_states: NUM_SOC_STATES,
            eta_c: 0.90,
            eta_d: 0.90,
            max_charge_mw: 2.0,
            max_discharge_mw: 1.8,
            degradation_cost: 2.0,
            terminal_ref_price: 45.0,
        }
    }

    #[inline]
    pub fn soc_mwh(soc_idx: usize) -> f64 {
        soc_idx as f64 * 0.1
    }

    #[inline]
    pub fn next_soc_idx(&self, soc_idx: usize, act: BatteryAction) -> Option<usize> {
        let next = soc_idx as isize + act.delta_soc_idx;
        if next >= 0 && (next as usize) < self.num_soc_states {
            Some(next as usize)
        } else {
            None
        }
    }

    #[inline]
    pub fn is_feasible(&self, soc_idx: usize, act: BatteryAction) -> bool {
        self.next_soc_idx(soc_idx, act).is_some()
    }

    #[inline]
    pub fn hourly_cashflow(&self, price: f64, act: BatteryAction) -> f64 {
        if act.power_mw > 0.0 {
            -price * act.power_mw - self.degradation_cost * act.power_mw
        } else if act.power_mw < 0.0 {
            let d = -act.power_mw;
            price * d - self.degradation_cost * d
        } else {
            0.0
        }
    }

    #[inline]
    pub fn terminal_value(&self, soc_idx: usize) -> f64 {
        self.terminal_ref_price * self.eta_d * Self::soc_mwh(soc_idx)
    }

    #[inline]
    pub fn equivalent_full_cycles(&self, stored_throughput_mwh: f64) -> f64 {
        stored_throughput_mwh / (2.0 * self.capacity_mwh)
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn test_battery_bounds_and_cashflow() {
        let bat = BatterySpec::default_spec();
        assert!((bat.terminal_value(bat.initial_soc_idx) - 81.0).abs() < 1e-12);
        assert!(bat
            .next_soc_idx(0, ACTIONS[ACTION_DISCHARGE_09MW])
            .is_none());
        assert!(bat.next_soc_idx(40, ACTIONS[ACTION_CHARGE_1MW]).is_none());
        assert!((bat.hourly_cashflow(-30.0, ACTIONS[ACTION_CHARGE_2MW]) - 56.0).abs() < 1e-12);
        assert!((bat.hourly_cashflow(200.0, ACTIONS[ACTION_DISCHARGE_18MW]) - 356.4).abs() < 1e-12);
    }
}
