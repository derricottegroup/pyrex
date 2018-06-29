def num_first_derivs(irc_energies,step_size):
    reaction_force = []
    for i in range(2,len(irc_energies)-2):
        current_force = -1.0*(-1.0*irc_energies[i+2] + 8.0*irc_energies[i+1] - 8.0*irc_energies[i-1] + 1.0*irc_energies[i-2])/(12.0*step_size)
        reaction_force.append(current_force)
    return reaction_force

def num_integrate(force_coordinates, reaction_force_values, lower_limit, upper_limit):
    integral_value = 0.0
    for i in range(lower_limit,upper_limit):
        integral_value += (force_coordinates[i+1] - force_coordinates[i])*((reaction_force_values[i] + reaction_force_values[i+1])/2.0)
    return integral_value
