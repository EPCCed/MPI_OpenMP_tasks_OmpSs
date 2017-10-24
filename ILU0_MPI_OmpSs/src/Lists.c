#include <stdio.h>
#include <stdlib.h>
#include "Lists.h"

// Verify if the List is empty or not
int emptyList (List lst) {
	return (lst.head == NULL);
}

// Push the data in a new node in the last position of lst
void PushList (ptr_List lst, void *data) {
	ptr_NodeList ptr = lst->end, nptr;

	if (data != NULL) {
		// Create the new node and fill its fields
		nptr = (ptr_NodeList) malloc (sizeof(NodeList));
		nptr->data = data;
		nptr->next = NULL;
		// Adjust the head and end fields of lst
		if (ptr == NULL)
			lst->head = nptr;
		else
			ptr->next = nptr;
		lst->end  = nptr;
	}
}

// Return the data included in the first position of lst and
// remove this node
void *PopList (ptr_List lst) {
	void *data = NULL;
	ptr_NodeList ptr = lst->head;

	if (ptr != NULL) {
		// Move the data before removing the node
		data = ptr->data;
		// Adjust the head and end fields of lst
		lst->head = ptr->next;
		if (lst->head == NULL)
			lst->end  = NULL;
		free (ptr);
	}

	return data;
}

// Remove the nodes of lst, if the application of tst on the 
// corresponding data returns true
void CleanList (ptr_List lst, TestNode tst) {
	ptr_NodeList ptr = lst->head, pptr = NULL;

	while (ptr != NULL) {
		// Ask if the node have to be removed
		if (tst (ptr->data)) {
			// Remove the node and move to the next one
			if (pptr == NULL) {
				lst->head = ptr->next;
				free (ptr);
				ptr = lst->head;
			} else {
				pptr->next = ptr->next;
				free (ptr);
				ptr = pptr->next;
			}
		}	else {
			// Move to the next node
			pptr = ptr;
			ptr  = pptr->next;
		}
	}
	// Adjust the head and end fields of lst
	if (pptr == NULL)
		lst->head  = pptr;
	lst->end  = pptr;
}

