function R_body_to_world = rpyToRotMatTest(rpy)
R_body_to_world = rx(rpy(1))' * ry(rpy(2))' * rz(rpy(3))';

