function z=eq_psi_y(x,y,C)
z = psi0_y(x,y,C(7))+C(1)*psi1_y(x,y,C(7))...
    +C(2)*psi2_y(x,y,C(6),C(7))+C(3)*psi3_y(x,y,C(6),C(7))...
    +C(4)*psi4_y(x,y,C(7))+C(5)*psi5_y(x,y,C(7))...
    +C(8)*psi6_y(x,y,C(7))+C(9)*psi7_y(x,y,C(7))...
    +C(10)*psi8_y(x,y,C(6),C(7))+C(11)*psi9_y(x,y,C(6),C(7))...
    +C(12)*psi10_y(x,y,C(7));