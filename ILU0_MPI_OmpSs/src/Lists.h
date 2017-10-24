#ifndef ListType

#define ListType 1

//Definition of the types related to the list datatype
typedef struct NodeList NodeList, *ptr_NodeList;
struct NodeList {
	void *data;
	ptr_NodeList next;
};
typedef struct List {
	NodeList *head, *end;
} List, *ptr_List;

// Verify if the List is empty or not
extern int emptyList (List lst);

// Push the data in a new node in the last position of lst
extern void PushList (ptr_List lst, void *data);

// Return the data included in the first position of lst and
// remove this node
extern void *PopList (ptr_List lst);

// Definition of the class of functions to clean a list
typedef int (*TestNode) (void *);

// Remove the nodes of lst, if the application of tst on the 
// corresponding data returns true
extern void CleanList (ptr_List lst, TestNode tst);

#endif
