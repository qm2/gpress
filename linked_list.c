#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include "linked_list.h"


int init(node **head, int column, int value)
{
    *head = malloc(sizeof(node));
    if (!*head) {
        fprintf(stderr, "Failed to init a linked list\n");
        return 1;
    }

    (*head)->column = column;
    (*head)->value = value;
    (*head)->next = NULL;
    return 0;
}

int insert(node **head, int column, int value)
{
    node *current = *head;
    // node *tmp;

    // do {
    //     tmp = current;
    //     current = current->next;
    // } while (current);

    /* create a new node after tmp */
    node *new = malloc(sizeof(node));
    if (!new) {
        fprintf(stderr, "Failed to insert a new element\n");
        return 1;
    }
    // new->next = NULL;
    // new->column = column;
    // new->value = value;

    // tmp->next = new;

    new->next = *head;
    new->column = column;
    new->value = value;
    *head= new;

    return 0;

}

void reverse(node **head)
{
    node *current = *head; 
    node *newnext = NULL;
    node *tmp;

    do {
        tmp = current->next;
        current->next = newnext;
        newnext = current;
        current = tmp;
    } while (current);

    *head = newnext;
    return;
}
