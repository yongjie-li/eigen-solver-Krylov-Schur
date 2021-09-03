function v_mapped = func_multiply_shifted_ABCD(lu_shift_A, lu_eqn1_coef, mat_B, mat_C, v_to_map)
    
    b_eqn1 = mat_C * klu(lu_shift_A, '\', v_to_map);
    x_eqn1 = klu(lu_eqn1_coef, '\', b_eqn1);
    v_mapped = klu(lu_shift_A, '\', v_to_map - mat_B * x_eqn1);
    
end