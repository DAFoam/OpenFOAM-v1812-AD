/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
     \\/     M anipulation  | Copyright (C) 2017 OpenCFD Ltd.
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

\*---------------------------------------------------------------------------*/

#include "error.H"
#include "nullObject.H"

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

template<class IteratorType>
inline const IteratorType& Foam::DLListBase::iterator_end()
{
    return *reinterpret_cast<const IteratorType*>(nullObjectPtr);
}


template<class IteratorType>
inline const IteratorType& Foam::DLListBase::iterator_rend()
{
    return *reinterpret_cast<const IteratorType*>(nullObjectPtr);
}


template<class IteratorType>
inline IteratorType Foam::DLListBase::iterator_first() const
{
    DLListBase* list = const_cast<DLListBase*>(this);

    if (size())
    {
        return IteratorType(list, const_cast<DLListBase::link*>(first_));
    }

    // Return an end iterator
    return IteratorType(list, nullptr);
}


template<class IteratorType>
inline IteratorType Foam::DLListBase::iterator_last() const
{
    DLListBase* list = const_cast<DLListBase*>(this);

    if (size())
    {
        return IteratorType(list, const_cast<DLListBase::link*>(last_));
    }

    // Return an end iterator
    return IteratorType(list, nullptr);
}


// * * * * * * * * * * * * * * * Iterator ends * * * * * * * * * * * * * * * //

inline const Foam::DLListBase::iterator& Foam::DLListBase::end()
{
    return iterator_end<DLListBase::iterator>();
}


inline const Foam::DLListBase::const_iterator&
Foam::DLListBase::cend() const
{
    return iterator_end<DLListBase::const_iterator>();
}


inline const Foam::DLListBase::const_iterator&
Foam::DLListBase::crend() const
{
    return iterator_rend<DLListBase::const_iterator>();
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

inline bool Foam::DLListBase::link::registered() const
{
    return prev_ != nullptr && next_ != nullptr;
}


inline void Foam::DLListBase::link::deregister()
{
    prev_ = next_ = nullptr;
}


inline Foam::label Foam::DLListBase::size() const
{
    return size_;
}


inline bool Foam::DLListBase::empty() const
{
    return !size_;
}


inline Foam::DLListBase::link*
Foam::DLListBase::first()
{
    if (!size_)
    {
        FatalErrorInFunction
            << "list is empty"
            << abort(FatalError);
    }
    return first_;
}


inline const Foam::DLListBase::link*
Foam::DLListBase::first() const
{
    if (!size_)
    {
        FatalErrorInFunction
            << "list is empty"
            << abort(FatalError);
    }
    return first_;
}


inline Foam::DLListBase::link*
Foam::DLListBase::last()
{
    if (!size_)
    {
        FatalErrorInFunction
            << "list is empty"
            << abort(FatalError);
    }
    return last_;
}


inline const Foam::DLListBase::link*
Foam::DLListBase::last() const
{
    if (!size_)
    {
        FatalErrorInFunction
            << "list is empty"
            << abort(FatalError);
    }
    return last_;
}


inline void Foam::DLListBase::clear()
{
    first_ = nullptr;
    last_  = nullptr;
    size_  = 0;
}


inline void Foam::DLListBase::swap(DLListBase& lst)
{
    std::swap(first_, lst.first_);
    std::swap(last_, lst.last_);
    std::swap(size_, lst.size_);
}


inline void Foam::DLListBase::transfer(DLListBase& lst)
{
    first_ = lst.first_;
    last_  = lst.last_;
    size_  = lst.size_;

    lst.clear();
}


inline Foam::DLListBase::link*
Foam::DLListBase::remove
(
    DLListBase::iterator& iter
)
{
    return remove(iter.node_);
}


inline Foam::DLListBase::link*
Foam::DLListBase::replace
(
    DLListBase::iterator& oldIter,
    DLListBase::link* newItem
)
{
    return replace(oldIter.node_, newItem);
}


// * * * * * * * * * * * * * * * STL iterator  * * * * * * * * * * * * * * * //

inline Foam::DLListBase::iterator::iterator
(
    DLListBase* list,
    DLListBase::link* item
)
:
    node_(item),
    list_(list),
    copy_()
{
    if (node_ != nullptr)
    {
        copy_ = *node_;
    }
}


inline Foam::DLListBase::link*
Foam::DLListBase::iterator::get_node() const
{
    return node_;
}


inline bool Foam::DLListBase::iterator::found() const
{
    return (node_ != nullptr);
}


inline void Foam::DLListBase::iterator::prev()
{
    if (list_)
    {
        // Check if the node_ is the first element (points to itself)
        // or if the list is empty because last element was removed
        if (node_ == copy_.prev_ || list_->first_ == nullptr)
        {
            node_ = nullptr;
        }
        else
        {
            node_ = copy_.prev_;
            copy_ = *node_;
        }
    }
}


inline void Foam::DLListBase::iterator::next()
{
    if (list_)
    {
        // Check if the node_ is the last element (points to itself)
        // or if the list is empty because last element was removed
        if (node_ == copy_.next_ || list_->last_ == nullptr)
        {
            node_ = nullptr;
        }
        else
        {
            node_ = copy_.next_;
            copy_ = *node_;
        }
    }
}


inline void Foam::DLListBase::iterator::operator=(const iterator& iter)
{
    node_ = iter.node_;
    list_ = iter.list_;
    copy_ = iter.copy_;
}


inline bool Foam::DLListBase::iterator::operator==(const iterator& iter) const
{
    return node_ == iter.node_;
}


inline bool Foam::DLListBase::iterator::operator!=(const iterator& iter) const
{
    return node_ != iter.node_;
}


inline Foam::DLListBase::iterator
Foam::DLListBase::begin()
{
    if (size())
    {
        return iterator_first<iterator>();
    }

    return end();
}


// * * * * * * * * * * * * * * STL const_iterator  * * * * * * * * * * * * * //

inline Foam::DLListBase::const_iterator::const_iterator
(
    const DLListBase* list,
    const DLListBase::link* item
)
:
    node_(item),
    list_(list)
{}


inline Foam::DLListBase::const_iterator::const_iterator
(
    const DLListBase::iterator& iter
)
:
    node_(iter.node_),
    list_(iter.list_)
{}


inline const Foam::DLListBase::link*
Foam::DLListBase::const_iterator::get_node() const
{
    return node_;
}


inline bool Foam::DLListBase::const_iterator::found() const
{
    return (node_ != nullptr);
}


inline void Foam::DLListBase::const_iterator::prev()
{
    if (list_ && node_)
    {
        if (node_ == list_->first_)
        {
            node_ = nullptr;
        }
        else
        {
            node_ = node_->prev_;
        }
    }
}


inline void Foam::DLListBase::const_iterator::next()
{
    if (list_ && node_)
    {
        if (node_ == list_->last_)
        {
            node_ = nullptr;
        }
        else
        {
            node_ = node_->next_;
        }
    }
}


inline bool Foam::DLListBase::const_iterator::operator==
(
    const const_iterator& iter
) const
{
    return node_ == iter.node_;
}


inline bool Foam::DLListBase::const_iterator::operator!=
(
    const const_iterator& iter
) const
{
    return node_ != iter.node_;
}


inline Foam::DLListBase::const_iterator
Foam::DLListBase::cbegin() const
{
    if (size())
    {
        return iterator_first<const_iterator>();
    }

    return cend();
}


inline Foam::DLListBase::const_iterator
Foam::DLListBase::crbegin() const
{
    if (size())
    {
        return iterator_last<const_iterator>();
    }

    return crend();
}


// ************************************************************************* //
