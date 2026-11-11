"""Read-only access to the existing docking leaderboard."""
from typing import Optional

from fastapi import APIRouter

from ..leaderboard import get_leaderboard

router = APIRouter()


@router.get("/api/leaderboard")
async def leaderboard(target_id: Optional[str] = None, top_n: int = 10):
    return get_leaderboard(target_id=target_id, top_n=top_n)
