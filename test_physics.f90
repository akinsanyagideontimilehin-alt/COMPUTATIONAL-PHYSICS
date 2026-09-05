program test_physics
    implicit none
    real :: position, velocity, time_step
    
    position = 0.0
    velocity = 5.0
    time_step = 0.5
    
    ! Simulate a basic linear update step
    position = position + velocity * time_step
    
    print *, "Testing compiler... Simulation step successful!"
    print *, "New Position:", position
end program test_physics
