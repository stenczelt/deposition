module neighbors

    implicit none

contains

    recursive subroutine find_neighbors(i, j, n_atoms, bonded, atom_visited, atom_belongs_to_cluster, cluster)

        implicit none

        integer :: i, j, k, n_atoms, cluster, atom_belongs_to_cluster(:), cluster_size
        logical :: atom_visited(:), bonded(:, :)

        if(.not. atom_visited(j) .and. i /= j .and. bonded(i, j))then
            atom_visited(j) = .true.
            atom_belongs_to_cluster(j) = cluster
            do k = 1, n_atoms
                call find_neighbors(j, k, n_atoms, bonded, atom_visited, atom_belongs_to_cluster, cluster)
            end do
        end if
    end subroutine


    subroutine get_distance(pos1, pos2, l1, l2, l3, d)

        !   Returns distance d between atoms i and j according
        !   to minimum image convention

        implicit none

        real*8 :: d, l1(1:2), l2(1:2), l3(1:2)
        real*8 :: x, y, z, lx, ly, lz, pos1(1:3), pos2(1:3)

        lx = l1(2) - l1(1)
        ly = l2(2) - l2(1)
        lz = l3(2) - l3(1)

        x = pos2(1) - pos1(1)
        y = pos2(2) - pos1(2)
        z = pos2(3) - pos1(3)

        if(dabs(x) > lx / 2.d0)then
            x = x - sign(lx, x)
        end if
        if(dabs(y) > ly / 2.d0)then
            y = y - sign(ly, y)
        end if
        if(dabs(z) > lz / 2.d0)then
            z = z - sign(lz, z)
        end if

        d = sqrt(x**2 + y**2 + z**2)

    end subroutine

end module
