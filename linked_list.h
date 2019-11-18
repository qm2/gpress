#ifndef LINKEDLIST_H
#define LINKEDLIST_H

//linked list element structure
typedef struct LinkedList {
	struct LinkedList* next; // Next element in case of a collision
    int column;
    int value;
}node;

int init(node **head, int column, int value);

int insert(node **head, int column, int value);

void reverse(node **head);

#endif
