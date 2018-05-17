MODULE basic_adjacency
USE fail
USE util
IMPLICIT NONE

TYPE :: adjacency_list
INTEGER, DIMENSION(:), ALLOCATABLE :: v
END TYPE

PUBLIC :: ALLOCATED
PUBLIC :: add_node
PUBLIC :: allocate_adj_list
PUBLIC :: allocate_array
PUBLIC :: ASSIGNMENT ( = )
PUBLIC :: deallocate_adj_list
PUBLIC :: deallocate_array
PUBLIC :: remove_node
PUBLIC :: SIZE

INTERFACE ALLOCATED
    MODULE PROCEDURE list_allocated
END INTERFACE

INTERFACE allocate_array
    MODULE PROCEDURE allocate_array_of_adj_lists
END INTERFACE

INTERFACE ASSIGNMENT ( = )
    MODULE PROCEDURE array_to_list
    MODULE PROCEDURE list_to_array
    MODULE PROCEDURE copy_list
END INTERFACE

INTERFACE deallocate_array
    MODULE PROCEDURE deallocate_array_of_adj_lists
END INTERFACE

INTERFACE SIZE
    MODULE PROCEDURE list_size
END INTERFACE

CONTAINS
    FUNCTION list_allocated ( list_1 )
        LOGICAL                          :: list_allocated
        TYPE(adjacency_list), INTENT(IN) :: list_1
        list_allocated = ALLOCATED ( list_1%v )
    END FUNCTION

    SUBROUTINE array_to_list ( list_1, arr )
        TYPE(adjacency_list),  INTENT(INOUT) :: list_1
        INTEGER, DIMENSION(:), INTENT(IN)    :: arr

        IF ( ALLOCATED ( list_1 ) ) THEN
            IF ( SIZE ( list_1 ) == SIZE ( arr ) ) THEN
                IF ( SIZE ( arr ) > 0 ) list_1%v = arr
            ELSE
                CALL deallocate_array ( list_1%v )
                CALL allocate_array ( list_1%v, SIZE ( arr ) )
                IF ( SIZE ( arr ) > 0 ) list_1%v = arr
            ENDIF
        ELSE
            CALL allocate_array ( list_1%v, SIZE ( arr ) )
            IF ( SIZE ( arr ) > 0 ) list_1%v = arr
        ENDIF
    END SUBROUTINE array_to_list

    SUBROUTINE list_to_array ( arr, list_1 )
        INTEGER, DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: arr
        TYPE(adjacency_list),               INTENT(IN)    :: list_1

        IF ( .NOT. ALLOCATED ( list_1 ) ) THEN
            CALL die ( 'Programming error: list_1 has not been allocated.', 'list_to_array' )
        ENDIF

        IF ( ALLOCATED ( arr ) ) THEN
            IF ( SIZE ( arr ) == SIZE ( list_1 ) ) THEN
                IF ( SIZE ( list_1 ) > 0 ) arr = list_1%v
            ELSE
                CALL deallocate_array ( arr )
                CALL allocate_array ( arr, SIZE ( list_1 ) )
                IF ( SIZE ( list_1 ) > 0 ) arr = list_1%v 
            ENDIF
       ELSE
           CALL allocate_array ( arr, SIZE ( list_1 ) )
           IF ( SIZE ( list_1 ) > 0 ) arr = list_1%v
       ENDIF
    END SUBROUTINE list_to_array

    SUBROUTINE copy_list ( list_1, list_2 )
        TYPE(adjacency_list), INTENT(INOUT) :: list_1
        TYPE(adjacency_list), INTENT(IN)    :: list_2

        IF ( .NOT. ALLOCATED ( list_2 ) ) THEN
            CALL die ( 'Programming error: list_2 has not been initialized.', 'copy_list' )
        ENDIF

        IF ( ALLOCATED ( list_1 ) ) THEN
            IF ( SIZE ( list_1 ) == SIZE ( list_2 ) ) THEN
                IF ( SIZE ( list_2 ) > 0 ) list_1%v = list_2%v
            ELSE
                CALL deallocate_array ( list_1%v )
                CALL allocate_array ( list_1%v, SIZE ( list_2%v ) )
                IF ( SIZE ( list_2 ) > 0 ) list_1%v = list_2%v
            ENDIF
        ELSE
            CALL allocate_array ( list_1%v, SIZE ( list_2%v ) )
            IF ( SIZE ( list_2 ) > 0 ) list_1%v = list_2%v
        ENDIF

    END SUBROUTINE copy_list

    FUNCTION list_size ( list_1 )
        INTEGER                          :: list_size
        TYPE(adjacency_list), INTENT(IN) :: list_1

        IF ( ALLOCATED ( list_1%v ) ) THEN
            list_size = SIZE ( list_1%v )
        ELSE
            CALL die ( 'Programming error: list_1 has not been allocated.', 'list_size' )
        ENDIF

    END FUNCTION list_size

    SUBROUTINE remove_node ( list_1, node_num )
        TYPE(adjacency_list), INTENT(INOUT) :: list_1
        INTEGER,              INTENT(IN)    :: node_num
!       Local vars:
        INTEGER, DIMENSION(:), ALLOCATABLE  :: vec
        INTEGER                             :: n 
        
!       Chekz:
        IF ( .NOT. ALLOCATED ( list_1 ) ) THEN
            CALL die ( 'Programming error: list_1 has not been allocated.', 'remove_node' )
        ENDIF

        IF ( SIZE ( list_1 ) == 0 ) THEN
            CALL warn ( 'Possible programming error. List has zero size. Cannot remove.', 'remove_node' )
            RETURN
        ENDIF

        n = COUNT ( list_1%v /= node_num )

!       Allocate temporary array:
        CALL allocate_array ( vec, n )

!       No change in vec size - no node found:
        IF ( n == SIZE ( list_1 ) ) THEN 
            WRITE(*,*) ' node to remove ', node_num, ' list ', list_1%v
            CALL die ( 'Programming error: no such node.', 'remove_node')
        ENDIF

        vec = PACK ( list_1%v, list_1%v /= node_num )
        list_1 = vec

!       Release temporary storage:
        CALL deallocate_array ( vec )

    END SUBROUTINE remove_node

    SUBROUTINE add_node ( list_1, node_num )
        TYPE(adjacency_list), INTENT(INOUT) :: list_1
        INTEGER,              INTENT(IN)    :: node_num
!       Local vars:
        INTEGER, DIMENSION(:), ALLOCATABLE  :: vec

        IF ( .NOT. ALLOCATED ( list_1 ) ) THEN
            CALL die ( 'Programming error. Cannot add node to unallocated list.', 'add_node' )
        ENDIF    

        CALL allocate_array ( vec, SIZE ( list_1) + 1 )
        vec(1:SIZE(list_1)) = list_1%v
        vec(SIZE(list_1)+1) = node_num
        CALL deallocate_array ( list_1%v )
        CALL allocate_array ( list_1%v, SIZE ( vec ) )
        list_1%v = vec
        CALL deallocate_array ( vec )

    END SUBROUTINE add_node
    
    SUBROUTINE allocate_adj_list ( list_1, list_size )
        TYPE(adjacency_list), INTENT(INOUT) :: list_1
        INTEGER,              INTENT(IN)    :: list_size

        IF ( list_size < 0 ) THEN
            CALL die ('Programming error: cannot allocate list_1 having negative length.', 'allocate_adj_list')
        ENDIF

        IF ( ALLOCATED ( list_1 ) ) THEN
            IF ( SIZE ( list_1 ) == list_size ) THEN
                RETURN
            ELSE
                CALL deallocate_array ( list_1%v )
            ENDIF
        ELSE
            CALL allocate_array ( list_1%v, list_size )
        ENDIF                  

        CALL allocate_array ( list_1%v, list_size )

   END SUBROUTINE allocate_adj_list

   SUBROUTINE deallocate_adj_list ( list_1 )
       TYPE(adjacency_list), INTENT(INOUT) :: list_1

       IF ( ALLOCATED ( list_1 ) ) CALL deallocate_array ( list_1%v )

   END SUBROUTINE deallocate_adj_list

   SUBROUTINE allocate_array_of_adj_lists ( arr, arr_size )
       TYPE(adjacency_list), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: arr
       INTEGER,                                         INTENT(IN)    :: arr_size
!      Local vars:
       INTEGER                                                        :: istat
!      Counters:
       INTEGER                                                        :: i

       ALLOCATE ( arr(arr_size), STAT=istat )
       IF ( istat == 0 ) THEN
           DO i = 1, arr_size

!              Since size of each list is unknown the only thing to do is to allocate them having zero length:
               CALL allocate_adj_list ( arr(i), 0 )
           ENDDO
       ELSE
           CALL die ( 'Failed to allocate array of adjacency lists.', 'allocate_array_of_adj_lists')
       ENDIF

    END SUBROUTINE allocate_array_of_adj_lists

    SUBROUTINE deallocate_array_of_adj_lists ( arr )
        TYPE(adjacency_list), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: arr
!       Local vars:
        INTEGER                                                        :: istat
!       Counters:
        INTEGER                                                        :: i

         IF ( ALLOCATED ( arr ) ) THEN

!            Deallocate everything in reverse order to allocation:
             DO i = 1, SIZE ( arr )
                 IF ( ALLOCATED ( arr(i) ) ) CALL deallocate_adj_list ( arr(i) )
             ENDDO


!            Final deallocation:
             DEALLOCATE ( arr, STAT=istat )
             IF ( istat /= 0 ) THEN
                 CALL die ('Failed to deallocate array of adjacency lists.', 'deallocate_array_of_adj_lists' )
             ENDIF
         ENDIF

    END SUBROUTINE deallocate_array_of_adj_lists 

    SUBROUTINE remove_node_from_adj_list ( arr, node )
        TYPE(adjacency_list), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: arr
        INTEGER,                                         INTENT(IN)    :: node
!       Local vars:
        INTEGER, DIMENSION(:), ALLOCATABLE                             :: list_of_neighbors
!       Counters:
        INTEGER                                                        :: i

!       Get list of neigbors for the node:
        list_of_neighbors = arr ( node )
        IF ( SIZE ( list_of_neighbors ) == 0 ) RETURN ! Nothing to remove

!       Remove node from lists of its neigbors:
        DO i = 1, SIZE ( list_of_neighbors )
            CALL remove_node (arr(list_of_neighbors(i)), node )
        ENDDO

!       Remove neighbors of node itself:
        DO i = 1, SIZE ( list_of_neighbors )
            CALL remove_node (arr(node), list_of_neighbors(i))
        ENDDO

!       Deallocate temporary array(s):
        CALL deallocate_array ( list_of_neighbors )

    END SUBROUTINE remove_node_from_adj_list

    SUBROUTINE add_node_to_adj_list ( arr, node, list_of_neighbors )
        TYPE(adjacency_list), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: arr
        INTEGER,                                         INTENT(IN)    :: node
        INTEGER,              DIMENSION(:), ALLOCATABLE, INTENT(IN)    :: list_of_neighbors
!       Local vars:
        INTEGER                                                        :: i

!       Recreate node itself:
        arr(node) = list_of_neighbors

!       Include node in lists of neighbors:      
        DO i = 1, SIZE ( list_of_neighbors )

!           Note that resulting adjacency list now unsorted:
            CALL add_node ( arr(list_of_neighbors(i)), node )
        ENDDO

    END SUBROUTINE add_node_to_adj_list
 
    SUBROUTINE check_neighbors_degree ( arr, list_of_neighbors, maxcon, ifail )
        TYPE(adjacency_list), DIMENSION(:), ALLOCATABLE, INTENT(INOUT) :: arr
        INTEGER,              DIMENSION(:), ALLOCATABLE, INTENT(IN)    :: list_of_neighbors
        INTEGER,                                         INTENT(IN)    :: maxcon 
        INTEGER,                                         INTENT(INOUT) :: ifail
!       Local vars:
        INTEGER                                                        :: i

        ifail = 0
        DO i = 1, SIZE ( list_of_neighbors )
            IF ( SIZE ( arr(list_of_neighbors(i)) ) > maxcon  ) THEN
                ifail = -1
                RETURN
            ENDIF  
        ENDDO
   END SUBROUTINE  check_neighbors_degree              

END MODULE basic_adjacency
