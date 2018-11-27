classdef PriorityQueue < handle
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% A priority queue class for MATLAB
%
% Author: Masashi Shimbo
%
% Description:
%
% This is a simple min-priority queue, equipped with a table needed for quick
% 'decreaeKey' operation.  Each stored item is assumed to be represented by an
% integer ID in range 1..maxSize, where maxSize is the maximum number of items
% to be stored (which must be declared at the time of object construction).
%
% Written for implementing Brandes's shortest-path betweenness algorithm,
% which appeared as Algorithm 4.2 in the following book:
% 
%   François Fouss, Marco Saerens and Masashi Shimbo (2016).
%   "Algorithms and models for network data and link analysis". 
%   Cambridge University Press.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
  properties (Access = private)
    last % number of elements in the heap
    heap % binary heap of [item ID, key (= priority)]
    pos % map from an item ID to its heap index
  end

  methods
    % constructor
    function self = PriorityQueue(maxSize)
      self.heap = zeros(maxSize, 2);
      self.pos = zeros(maxSize, 1);
      self.last = 0;
    end

    function clear(self)
      self.heap = zeros(size(self.heap));
      self.pos = zeros(size(self.pos));
      self.last = 0;
    end

    function sz = size(self)
      sz = self.last;
    end

    % insert a new item
    function insert(self, item, key)
      if self.last == size(self.heap, 1)
	error('priority queue full');
      end

      if item < 0 || item > size(self.heap, 1)
	error('item out of range');
      end
      % place new item at the end
      self.last = self.last + 1;
      self.heap(self.last, :) = [item key];
      self.pos(item) = self.last;
      % heapify
      self.upHeap(self.last);
    end

    % return the item with the minimum key
    function [item, key] = min(self)
      if self.last == 0
	error('priority queue empty');
      end
      item = self.heap(1, 1);
      key  = self.heap(1, 2);
    end
    
    % extract the item at the top of the heap (with the minimum key value)
    % the item is removed from the heap
    function [item, key] = extractMin(self)
      [item, key] = min(self);
      
      % clear records for the removed item
      self.pos(item) = 0;

      % overwrite position 1 in the heap with the last item
      if self.last == 1
	self.last = 0;
      else
        % bring the last item to the vacant top
	self.heap(1, :) = self.heap(self.last, :);
	self.last = self.last - 1;
	self.pos(self.heap(1)) = 1;
        % heapify	
	self.downHeap(1);
      end
    end

    % update the key value of an item
    function decreaseKey(self, item, key)
      if item < 0 || item > size(self.heap, 1)
	error('item (%d) out of range', item)
      end
      j = self.pos(item);
      if j == 0
	error('item (%d) not found in priority queue', item);
      end
      oldKey = self.heap(self.pos(item), 2);
      if key >= oldKey
	error('new key (%f) not smaller than current key (%f)', key, oldKey);
      end
      self.heap(j, 2) = key;
      self.upHeap(j);
    end

    % return key associated with an item
    function k = key(self, item)
      j = self.pos(item);
      if j == 0
	error('item (%d) not found in priority queue', item);
      end
      k = self.heap(j, 2);
    end
  end

  methods (Access = private)
    function upHeap(self,  i)
      while i > 1
	j = floor(i / 2); % j = parent of i
	if self.heap(j, 2) <= self.heap(i, 2) % if already heapified
	  break
	end
	self.swap(i, j);
	i = j;
      end
    end

    function downHeap(self,  i)
      j = 2 * i;
      while j < self.last % strict inequality ensures that i has two children (j and j+1)
	if self.heap(j + 1, 2) < self.heap(j, 2)
	  j = j + 1; % j now points to the smallest of the left and right children
	end
	if self.heap(i, 2) <= self.heap(j, 2) % if already heapified
	  break
	end
	self.swap(i, j); % otherwise, push current key at i down (to j)
	i = j;
	j = 2 * i;
      end
      % corner case, when node i has a single child j (== last in the heap)
      if j == self.last && self.heap(j, 2) < self.heap(i, 2)
	self.swap(i, j);
      end
    end

    % exchange two heap entries
    function swap(self, i, j)
      self.heap([i j], :) = self.heap([j i], :);
      self.pos(self.heap([i j], 1)) = [i j];
    end
  end
end
