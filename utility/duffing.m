function qdot=duffing(t,q);
qdot=[q(2);
    -0.1*q(2)+q(1)-q(1)^3];