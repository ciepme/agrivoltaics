function npv = calc_discounted_cash_flow(income, discount_rate, investigation_period)
    exponent_value = 0:investigation_period-1;
    npv = sum(income ./ (1 + discount_rate).^exponent_value);
end