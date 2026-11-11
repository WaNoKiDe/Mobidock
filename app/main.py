"""MOBIDock FastAPI entry point."""
import logging
import os
from contextlib import asynccontextmanager

from fastapi import FastAPI
from fastapi.middleware.cors import CORSMiddleware
from fastapi.responses import FileResponse, JSONResponse
from fastapi.staticfiles import StaticFiles
from fastapi.routing import APIRoute
from starlette.routing import Mount

from .config import STATIC_DIR
from .routes import dock, drug_warnings, insight, leaderboard, ligands, pipeline, similarity, synth, targets, viewer

logging.basicConfig(level=logging.INFO, format="%(asctime)s  %(levelname)-8s %(name)s: %(message)s", datefmt="%H:%M:%S")
logger = logging.getLogger("mobidock")


@asynccontextmanager
async def lifespan(app: FastAPI):
    logger.info("MOBIDock is ready — http://localhost:8000")
    yield


_public_mode = os.environ.get("MOBIDOCK_PUBLIC_MODE", "").strip().lower() in {"1", "true", "yes", "on"}
if _public_mode:
    from gradio import Server
    app = Server()
else:
    app = FastAPI(title="MOBIDock", version="1.0.0", lifespan=lifespan)
app.add_middleware(CORSMiddleware, allow_origins=["*"], allow_methods=["*"], allow_headers=["*"])
app.mount("/static", StaticFiles(directory=str(STATIC_DIR)), name="static")


@app.middleware("http")
async def _no_cache(request, call_next):
    resp = await call_next(request)
    if request.url.path == "/" or request.url.path.startswith("/static/"):
        resp.headers["Cache-Control"] = "no-store"
    return resp


for _r in (dock.router, targets.router, ligands.router, leaderboard.router, similarity.router, synth.router, drug_warnings.router, viewer.router, insight.router, pipeline.router):
    if _public_mode:
        # Register each route through FastAPI itself. The recent FastAPI router
        # stores include_router() calls lazily, and Gradio's launch lifespan
        # does not always expand those deferred nodes before the first request.
        for _route in _r.routes:
            if isinstance(_route, APIRoute):
                app.add_api_route(
                    _route.path, _route.endpoint, methods=_route.methods,
                    response_model=_route.response_model, status_code=_route.status_code,
                    tags=_route.tags, dependencies=_route.dependencies, summary=_route.summary,
                    description=_route.description, response_description=_route.response_description,
                    responses=_route.responses, deprecated=_route.deprecated, operation_id=_route.operation_id,
                    response_model_include=_route.response_model_include,
                    response_model_exclude=_route.response_model_exclude,
                    response_model_by_alias=_route.response_model_by_alias,
                    response_model_exclude_unset=_route.response_model_exclude_unset,
                    response_model_exclude_defaults=_route.response_model_exclude_defaults,
                    response_model_exclude_none=_route.response_model_exclude_none,
                    include_in_schema=_route.include_in_schema, response_class=_route.response_class,
                    name=_route.name, openapi_extra=_route.openapi_extra,
                )
            elif isinstance(_route, Mount):
                app.mount(_route.path, _route.app, name=_route.name)
    else:
        app.include_router(_r)


async def _serve_index():
    index = STATIC_DIR / "index.html"
    if not index.exists():
        return JSONResponse({"error": "Frontend not found"}, status_code=404)
    return FileResponse(str(index))


if _public_mode:
    # Serve the custom MobiDock frontend directly. Gradio's own frontend is
    # disabled in the public launcher so it cannot catch these paths first.
    app.add_api_route("/", _serve_index, methods=["GET"], include_in_schema=False)
    app.add_api_route("/mobidock", _serve_index, methods=["GET"], include_in_schema=False)
else:
    app.add_api_route("/", _serve_index, methods=["GET"], include_in_schema=False)


from . import results_export
if _public_mode:
    for _route in results_export.router.routes:
        if isinstance(_route, APIRoute):
            app.add_api_route(
                _route.path, _route.endpoint, methods=_route.methods,
                response_model=_route.response_model, status_code=_route.status_code,
                tags=_route.tags, dependencies=_route.dependencies, summary=_route.summary,
                description=_route.description, response_description=_route.response_description,
                responses=_route.responses, deprecated=_route.deprecated, operation_id=_route.operation_id,
                response_model_include=_route.response_model_include,
                response_model_exclude=_route.response_model_exclude,
                response_model_by_alias=_route.response_model_by_alias,
                response_model_exclude_unset=_route.response_model_exclude_unset,
                response_model_exclude_defaults=_route.response_model_exclude_defaults,
                response_model_exclude_none=_route.response_model_exclude_none,
                include_in_schema=_route.include_in_schema, response_class=_route.response_class,
                name=_route.name, openapi_extra=_route.openapi_extra,
            )
        elif isinstance(_route, Mount):
            app.mount(_route.path, _route.app, name=_route.name)
else:
    app.include_router(results_export.router)
